#pragma once

#include "Scalar.h"
#include "VectorMath.h"

#include <fstream>

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

class SimplePolygonMesh {
public:
  SimplePolygonMesh();
  
  // == Mesh data
  std::vector<std::vector<size_t>> faces;
  std::vector<std::vector<Scalar>> vertices;
  std::vector<std::vector<std::vector<Scalar>>> uv; // optional UV coords, in correspondence with faces array

  // == Accessors
  inline size_t nFaces() const { return faces.size(); }
  inline size_t nVertices() const { return vertices.size(); }
  inline bool hasParameterization() const { return !uv.empty(); }

  // == Mutators

  // Empty all data arrays
  void clear();


  // === Input & ouput

  void readMeshFromFile(std::istream& in, std::string type);
  void readMeshFromFile(std::string filename, std::string type = "");
  void readMeshFromFile(std::string filename, std::string type,
                        std::string& detectedType); // also returns type string for filetype that was used
  void writeMesh(std::ostream& out, std::string type);
  void writeMesh(std::string filename, std::string type = "");

  void MinMaxBounds(std::vector<Scalar> &min_v, std::vector<Scalar> &max_v);
  void ScaleAndCenter();


private:
  std::string detectFileType(std::string filename);

  // Read helpers
  void readMeshFromObjFile(std::istream& in);

  // Write helpers
  void writeMeshObj(std::ostream& out);
};

} // namespace cpm
