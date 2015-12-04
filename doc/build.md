[OGDF](/README.md) » Build Guide

# Build Guide

The OGDF build configuration is generated using [CMake](http://www.cmake.org/).

## Requirements

 * CMake 3.1+
 * C++11 compliant compiler
   * gcc 4.8+
   * clang 3.5+
   * Microsoft Visual C++ 2015+
 * GNU Make (in most cases)
 * Doxygen (optional)

## Build Configuration

CMake is a meta-build system that will generate your actual build system.
The most common build systems include Unix Makefiles, Visual Studio project files, and XCode project files.
We refer to [Running CMake](http://www.cmake.org/runningcmake/) for extensive information on using CMake.

Note that CMake allows you to place the generated build system in a separate folder, thus allowing several parallel build configurations (called out-of-source builds). We recommend following this approach.

## Unix Examples

All of the examples cover unix systems. On windows, CMake provides a self-explanatory graphical user interface.

We assume that the directory `~/OGDF` contains a clone of this repository.
In each of the examples we are initially located in the home directory `~`.

CMake will present you with several options for configuring your build.
All of these options are explained by the CMake interface.

### Default Configuration

This will generate the default build system inside the OGDF directory.
Such a configuration is sufficient if you want to compile the Library with a single configuration only.

```
$ cd OGDF
$ cmake .
-- Configuring done
-- Generating done
-- Build files have been written to: ~/OGDF
```

Assuming your system supports Unix Makefiles you could now start the actual build:
Note that invoking make with the `-j` flag will execute jobs in parallel thus speeding up your build.

```
$ make -j8
Scanning dependencies of target COIN
[  0%] Building CXX object CMakeFiles/COIN.dir/src/coin/Symphony/cp_wrapper.cpp.o
[  0%] Building CXX object CMakeFiles/COIN.dir/src/coin/Symphony/lp.cpp.o
[  0%] Building CXX object CMakeFiles/COIN.dir/src/coin/Symphony/lp_wrapper.cpp.o
[  0%] Building CXX object CMakeFiles/COIN.dir/src/coin/Symphony/master_prep.cpp.o
[  0%] Building CXX object CMakeFiles/COIN.dir/src/coin/Symphony/cp_proccomm.cpp.o
[  1%] Building CXX object CMakeFiles/COIN.dir/src/coin/Symphony/master_prep_sr.cpp.o
...
```

### Out-of-Source Build

An out-of-source build does not modify the source directory.
Thus, generating a second out-of-source build will not create any conflicts with the previous one.

```
$ mkdir my-out-of-source-build
$ cd my-out-of-source-build
$ cmake ../OGDF
-- Configuring done
-- Generating done
-- Build files have been written to: ~/my-out-of-source-build
```

### Custom Configuration

Creating multiple build configurations is of no use if the configurations do not differ.
Run `ccmake` instead of `cmake` to modify the build configuration before generating the build system.

This includes the specification of the linear program solver, whether to include OGDF specific assertions, additional include directories and the likes.
Running `ccmake` also allows you to specify the compiler, linker and the respective flags.

```
$ cd OGDF
$ ccmake .

BUILD_SHARED_LIB              > OFF
CMAKE_BUILD_TYPE                 Release
CMAKE_INSTALL_PREFIX             /usr/local
COIN_EXTERNAL_SOLVER_INCLUDE_D
COIN_EXTERNAL_SOLVER_LIBRARIES
COIN_SOLVER                      CLP
OGDF_DEBUG                       ON
OGDF_MEMORY_MANAGER              POOL_TS
OGDF_SEPARATE_TESTS              OFF
CMAKE_CXX_COMPILER               /usr/bin/clang++
CMAKE_CXX_FLAGS_DEBUG            -g
CMAKE_CXX_FLAGS_RELEASE          -O3 -DNDEBUG

BUILD_SHARED_LIBS: Whether to build shared libraries instead of static ones.
Press [enter] to edit option
Press [c] to configure       Press [g] to generate and exit
Press [h] for help           Press [q] to quit without generating
Press [t] to toggle advanced mode (Currently Off)
```
