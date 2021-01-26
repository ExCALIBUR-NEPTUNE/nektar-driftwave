# Sample Hasegawa-Wakatani solver

This repo contains a sample solver for the Hasegawa-Wakatani equations based on
Nektar++.

More details on the formulation are [available in the brief documentation
workbook](doc.ipynb).


## Compilation & installation

The easiest way to build this solver is via Docker; from the top directory, build an image by running
```
docker build -t nektar-hw-solver .
```
and then run an interactive shell with
```
docker run -it nektar-hw-solver /bin/bash
```
If you have compiled and installed Nektar++ elsewhere, you can manually configure and compile this solver through CMake:
```
mkdir build && cd build
cmake -DNektar++_DIR=/path/to/lib/nektar++/cmake ..
make install
```
By default and unless a `CMAKE_INSTALL_PREFIX` is supplied binaries will appear in `build/dist/bin`. If Nektar++ is built with MPI, this will be supported in the `DriftWaveSolver`.
