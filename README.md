# Nektar++ solver for Hasegawa-Wakatani equations

This repo contains a sample solver for the Hasegawa-Wakatani equations using Nektar++.

More details on the formulation are [available in the brief documentation
workbook](doc.ipynb).


## Compilation & installation

The easiest way to build this solver is via Docker; from the top directory, build an image by running
```
docker build -t nektar-driftwave .
```
and then run an interactive shell with
```
docker run -it nektar-driftwave /bin/bash
```
If you have compiled and installed Nektar++ elsewhere, you can manually configure and compile this solver through CMake:
```
mkdir build && cd build
cmake -DNektar++_DIR=/path/to/lib/nektar++/cmake ..
make install
```
By default and unless a `CMAKE_INSTALL_PREFIX` is supplied binaries will appear in `build/dist/bin`. If Nektar++ is built with MPI, this will be supported in the `DriftWaveSolver`.

## Running the enclosed example

Assuming that you are working from the shell in the Docker image, you can run the enclosed example either in serial or parallel as follows:
```
cp /usr/local/share/doc/nektar++/nektar-driftwave/*.xml .
DriftWaveSolver driftwave.xml square_quads.xml               # For serial execution
mpirun -n 20 DriftWaveSolver driftwave.xml square_quads.xml  # For parallel execution
```
