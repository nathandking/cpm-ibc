C++ implementation of the closest point method for surface PDEs with interior (and exterior) boundary conditions. This code accompanies the paper "A Closest Point Method for PDEs on Manifolds with Interior Boundary Conditions for Geometry Processing".

---

## Building

There are some dependencies: [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page), [Polyscope](https://polyscope.run/), [LBFGS++](https://lbfgspp.statr.me/), [fcpw](https://github.com/rohan-sawhney/fcpw). You must download Eigen and cmake will use find_package() to include it. Polyscope, LBFGS++, and fcpw will be downloaded when you clone the repository using the following command:
```
git clone --recurse-submodules https://github.com/nathandking/cpm-ibc.git
```

ccmake is used to build the project. You should simply need to run the following commands: 
```
mkdir build 
cd build/ 
ccmake .. #choose the desired options (see below), configure, and generate
make 
```
---

## Running Examples

Once you build, the examples in the examples/ and convergence_studies/ directories will reside in build/bin/. To run an example simply run the name as a command. For example:
```
cd bin/
DiffusionCurvesCodimZero
```
---

## Options

BUILD_CUSTOM_SOLVER: If ON, use our partially matrix-free linear system solver. If OFF, option of BUILD_EIGEN_SPARSELU will appear.

BUILD_EIGEN_SPARSELU: If ON, use Eigen's direct solver, SparseLU. If OFF, use Eigen's BiCGSTAB iterative solver.

BUILD_ENABLE_SPARSE_GRID_SUPPORT: (Note: you must also set USE_SPARSE_GRID to ON). If ON, use our memory-efficient sparse-grid construction of the computational tube. If OFF, use our less memory-efficient, but faster construction of the computational tube.

USE_POLYSCOPE: If ON, visualization will occur using [Polyscope](https://polyscope.run/). If OFF, no visualization will occur (useful for running convergence_studies/ or on remote servers).

All other options are for included libraries.

---
Authors: [Nathan King](https://nathandking.github.io/) and [Haozhe Su](https://soldierdown.github.io/) (custom solver and sparse-grid support).

If cpm-interior-boundary-conditions contributes to an academic publication, cite it as:
```bib
@misc{cpm-ibc,
  title = {C++ Code for the Closest Point Method with Interior Boundary Conditions},
  author = {Nathan King and Haozhe Su},
  note = {https://github.com/nathandking/cpm-ibc},
  year = {2024}
}
```