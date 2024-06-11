#include <iterator>
#include <iomanip>
#include "SimplePolygonMesh.h"

namespace cpm {


// This code has been copied (and modified) from geometry-central, see licesne information below
// MIT License

// Copyright (c) 2017-2019 Nicholas Sharp and the geometry-central contributors

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.


SimplePolygonMesh::SimplePolygonMesh() {}

namespace { // helpers for parsing

class Index {
public:
  Index() {}

  Index(long long int v, long long int vt, long long int vn) : position(v), uv(vt), normal(vn) {}

  bool operator<(const Index& i) const {
    if (position < i.position) return true;
    if (position > i.position) return false;
    if (uv < i.uv) return true;
    if (uv > i.uv) return false;
    if (normal < i.normal) return true;
    if (normal > i.normal) return false;

    return false;
  }

  long long int position = -1;
  long long int uv = -1;
  long long int normal = -1;
};

Index parseFaceIndex(const std::string& token) {
  std::stringstream in(token);
  std::string indexString;
  int indices[3] = {1, 1, 1};

  int i = 0;
  while (std::getline(in, indexString, '/')) {
    if (indexString != "\\") {
      std::stringstream ss(indexString);
      ss >> indices[i++];
    }
  }

  // decrement since indices in OBJ files are 1-based
  return Index(indices[0] - 1, indices[1] - 1, indices[2] - 1);
}

std::vector<std::string> supportedMeshTypes = {"obj"};

} // namespace


std::string SimplePolygonMesh::detectFileType(std::string filename) {
  std::string::size_type sepInd = filename.rfind('.');
  std::string type;

  if (sepInd != std::string::npos) {
    std::string extension;
    extension = filename.substr(sepInd + 1);

    // Convert to all lowercase
    std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
    type = extension;
  } else {
    throw std::runtime_error("Could not auto-detect file type to load mesh from " + filename);
  }

  // Check if this is one of the filetypes we're aware of
  if (std::find(std::begin(supportedMeshTypes), std::end(supportedMeshTypes), type) == std::end(supportedMeshTypes)) {
    throw std::runtime_error("Detected file type " + type + " to load mesh from " + filename +
                             ". This is not a supported file type.");
  }

  return type;
}

void SimplePolygonMesh::readMeshFromFile(std::string filename, std::string type) {
  std::string unused;
  readMeshFromFile(filename, type, unused);
}

void SimplePolygonMesh::readMeshFromFile(std::string filename, std::string type, std::string& detectedType) {

  // Attempt to detect filename
  bool typeGiven = type != "";
  if (!typeGiven) {
    type = detectFileType(filename);
  }

  // == Open the file and load it
  // NOTE: Intentionally always open the stream as binary, even though some of the subsequent formats are plaintext and
  // others are binary.  The only real difference is that non-binary mode performs automatic translation of line ending
  // characters (e.g. \r\n --> \n from DOS). However, this behavior is platform-dependent and having platform-dependent
  // behavior seems more confusing then just handling the newlines properly in the parsers.
  std::ifstream inStream(filename, std::ios::binary);
  if (!inStream) throw std::runtime_error("couldn't open file " + filename);
  readMeshFromFile(inStream, type);

  detectedType = type;
}

void SimplePolygonMesh::readMeshFromFile(std::istream& in, std::string type) {

  if (type == "obj") {
    readMeshFromObjFile(in);
  } else {
    throw std::runtime_error("Did not recognize mesh file type " + type);
  }
}

// Read a .obj file containing a polygon mesh
void SimplePolygonMesh::readMeshFromObjFile(std::istream& in) {
  clear();

  // corner UV coords, unpacked below
  std::vector<std::vector<Scalar>> coords;
  std::vector<std::vector<size_t>> polygonCoordInds;

  // parse obj format
  std::string line;
  while (getline(in, line)) {
    std::stringstream ss(line);
    std::string token;

    ss >> token;

    if (token == "v") {
      std::vector<Scalar> position((std::istream_iterator<Scalar>(ss)), 
                    std::istream_iterator<Scalar>());

      vertices.push_back(position);

    } else if (token == "vt") {
      double u, v;
      ss >> u >> v;

      coords.push_back(std::vector<Scalar>{u, v});

    } else if (token == "vn") {
      // Do nothing

    } else if (token == "f") {
      std::vector<size_t> face;
      std::vector<size_t> faceCoordInds;
      while (ss >> token) {
        Index index = parseFaceIndex(token);
        if (index.position < 0) {
          getline(in, line);
          size_t i = line.find_first_not_of("\t\n\v\f\r ");
          index = parseFaceIndex(line.substr(i));
        }

        face.push_back(index.position);

        if (index.uv != -1) {
          faceCoordInds.push_back(index.uv);
        }
      }

      faces.push_back(face);
      if (!faceCoordInds.empty()) {
        polygonCoordInds.push_back(faceCoordInds);
      }
    } else if (token == "l") {
      std::vector<size_t> inds;
      std::istringstream lineStream(line.substr(2)); // all but first char and space

      // parse out list of indices
      std::string token;
      while (std::getline(lineStream, token, ' ')) {
        std::stringstream tokenStream(token);
        size_t i;
        tokenStream >> i;
        inds.push_back(i - 1);
      }

      for(size_t i = 0; i < inds.size() - 1; ++i)
      {
          std::vector<size_t> edge(2);
          edge[0] = inds[i];
          edge[1] = inds[i + 1];
          faces.push_back(edge);
      }
    }
  }

  // If we got uv coords, unpack them in to per-corner values
  if (!polygonCoordInds.empty()) {
    for (std::vector<size_t>& faceCoordInd : polygonCoordInds) {
      uv.emplace_back();
      std::vector<std::vector<Scalar>>& faceCoord = uv.back();
      for (size_t i : faceCoordInd) {
        if (i < coords.size()) faceCoord.push_back(coords[i]);
      }
    }
  }
}


void SimplePolygonMesh::clear() {
  faces.clear();
  vertices.clear();
  uv.clear();
}


void SimplePolygonMesh::writeMesh(std::string filename, std::string type) {

  // Auto-detect type if needed
  bool typeGiven = type != "";
  if (!typeGiven) {
    type = detectFileType(filename);
  }

  // NOTE if/when we ever start writing binary file formats, will need to open in binary mode for those formats
  std::ofstream outStream(filename);
  if (!outStream) throw std::runtime_error("couldn't open output file " + filename);
  writeMesh(outStream, type);
}

void SimplePolygonMesh::writeMesh(std::ostream& out, std::string type) {
  if (type == "obj") {
    return writeMeshObj(out);
  } else {
    throw std::runtime_error("Write mesh file type " + type + " not supported");
  }
}

void SimplePolygonMesh::writeMeshObj(std::ostream& out) {

  // Make sure we write out at full precision
  out << std::setprecision(std::numeric_limits<double>::max_digits10);

  // Write header
  out << "#  vertices: " << vertices.size() << std::endl;
  if(faces[0].size() == 3)
  {
    out << "#     faces: " << faces.size() << std::endl;
  }
  else if(faces[0].size() == 2)
  {
    out << "#     edges: " << faces.size() << std::endl;
  }
  out << std::endl;

  // Write vertices
  for (std::vector<Scalar> p : vertices) {
    out << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
  }


  // Write texture coords (if present)
  for (std::vector<std::vector<Scalar>>& coords : uv) {
    for (std::vector<Scalar> c : coords) {
      out << "vt " << c[0] << " " << c[1] << std::endl;
    }
  }

  // Write faces
  if(faces[0].size() == 3)
  {
    size_t iC = 0;
    for (std::vector<size_t>& face : faces) {
      out << "f";
      for (size_t ind : face) {
        out << " " << (ind + 1);

        if (!uv.empty()) {
          out << "/" << (iC + 1);
          iC++;
        }
      }
      out << std::endl;
    }
  }
  else if(faces[0].size() == 2)
  {
    out << "l";
    for (std::vector<size_t>& face : faces) {
      out << " " << (face[0] + 1);
    }
    out << " " << (faces[nFaces() - 1][1] + 1);
    out << std::endl;
  }
}


// the min and max of the mesh in each direction
void SimplePolygonMesh::MinMaxBounds(std::vector<Scalar> &min_v, std::vector<Scalar> &max_v)
{
  min_v.assign(3, std::numeric_limits<Scalar>::max());
  max_v.assign(3, std::numeric_limits<Scalar>::min());
  for(size_t i = 0; i < nVertices(); ++i)
  {
      for(size_t d = 0; d < 3; ++d)
      {
          if(vertices[i][d] < min_v[d])
          {
              min_v[d] = vertices[i][d];
          }

          if(vertices[i][d] > max_v[d])
          {
              max_v[d] = vertices[i][d];
          }
      }
  }
}



// scale mesh to be in [-1,1]^3 and centered about (0,0,0)
void SimplePolygonMesh::ScaleAndCenter()
{
  std::vector<Scalar> min_v, max_v;
  MinMaxBounds(min_v, max_v);

  std::vector<Scalar> diff = max_v - min_v;

  //find dimension that has the biggest diff
  std::vector<size_t> ind = SortedIndices(diff);
  size_t max_diff_ind = ind[2];

  // shift object to center about (0,0,0)
  for(size_t i = 0; i < nVertices(); ++i)
  {
    for(size_t d = 0; d < 3; ++d)
    {
      vertices[i][d] -= (min_v[d] + 0.5 * diff[d]);
    }
  }

  // scale object to be in [-1,1]^3
  for(size_t i = 0; i < nVertices(); ++i)
  {
    for(size_t d = 0; d < 3; ++d)
    {
      vertices[i][d] /= (0.5 * diff[max_diff_ind]);
    }
  }
}


} // namespace cpm

