[![CS440 Banner](https://rgl.s3.eu-central-1.amazonaws.com/media/uploads/wjakob/2017/02/16/cs440-logo_web.jpg)](https://rgl.s3.eu-central-1.amazonaws.com/media/uploads/wjakob/2017/02/20/cs440-rgl.jpg)

## LumenTree

LumenTree is a portable ray tracer written in C++. The implementation is based on [Nori](https://wjakob.github.io/nori), an simple educational ray tracer used in the course [Advanced Computer Graphics](https://rgl.epfl.ch/courses/ACG17) taught at EPFL. The program supports Windows, Linux, and Mac OS environment and provides basic functionality of simple ray tracing.

---

### Cloning and Compilation
#### Cloning
LumenTree uses submodules to manage external libraries. To clone from github, add `--recursive` to command line.
```
git clone --recursive https://github.com/li779/LumenTree.git
```
Cloning repository without `--recursive` will result in empty folders, if so, try run 
```
git submodule update --init --recursive
```
to get all external libraries
#### Compilation
Before compiling the project, make sure have CMake installed. On Linux/macOS, make sure have a gcc compiler; On Windows, a Visual Studio 2019 will work.
On Linux/macOS, create a building directory `build` and compile project into build
```
$ cd path-to-LumenTree\
$ mkdir build
$ cd build
$ cmake ..
$ make -j
``` 
On Windows, also create a building directory `build`. Use CMake to compile project into build. Make sure select "Visual Studio 16 2019" to generate visual studio project file. After that, open `LumenTree.sln` with Visual Studio and compile.

---

### Implemented Functionality
#### Accelation Data Structure
1. LumenTree natively uses a simple **octree** to store all shapes in the scene
2. A high performance Intel Embree ray tracing kernal can be used by changing CMake options

#### BSDF
1. diffusive material
2. dielectric material
3. perfect mirror material

#### Integrator
1. simple normal integrator
2. whitted integrator

#### Light
1. area light for any object in scene

---
### Known Issues
There is a known issue with the NanoGUI version that Nori uses: on Linux systems with an integrated Intel GPU, a bug in the Mesa graphics drivers causes the GUI to freeze on startup. A workaround is to temporarily switch to an older Mesa driver to run Nori. This can be done by running.
```
export MESA_LOADER_DRIVER_OVERRIDE=i965
```
