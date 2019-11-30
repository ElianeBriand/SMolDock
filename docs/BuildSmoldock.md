
The build environement for SmolDock is  honestly somewhat cumbersome, because you may need to build a lot of dependencies 
from sources in some cases (eg wanting fully static executables). 

Further more there is manual modification of the CMakeList.txt file.

While we hope to simplify this process in the future, do know that once the work is done you will be able to build new versions
with only incremental changes, if you keep a copy of the "build setup section" (see below) to paste in future release.

Building on linux is supported. It should build on other UNIX-like, possibly with minor modification to the modus operandi.
 It likely will not build on Windows (if only due to the libunwind dependency). 

# Fetching and building the dependencies

Some of SmolDock's dependencies are shipped with the sources (due to slight modifications). Other you must install yourself :
- RDKit (BSD 3-clauses)
- Vc (BSD 3-clauses)
- libunwind (MIT licence)
- Eigen (MPL2 licence)
- Python libraries for the python modules (PSF)
- Boost (Boost software licence)
- Any standard-compliant MPI implementation (tested on OpenMPI)

For all but RDKit, you may want to use your distribution's package (often the version with "-dev" appended), or official build. However, a specific build process
is necessary for RDKit, and may be necessary for other libraries if they aren't available/suitable in your distribution package
manager repository.

On the other hand, you may have specific requirement for the resulting SmolDock build, like that it be fully static. For that,
you may have to do additional work detailled in the next sections.

In addition, a particular modification of Intel TBB library, namely the removal of a lock in initialization, may cause platform
or executable specific problem. If that's the case, you may want to follow the dedicated sections instruction here and afterward 
("If needed : External Intel TBB")

Make sure that you install the package wich containes the headers and static library for all of these dependencies. (often the 
package is split in "libfoo" and "libfoo-dev", and you need both to be able to build)

## Eigen
**General case** : install through your distribution package manager to the system-appropriate directory (/usr, typically)

**Special case** : not applicable (?)

## Python (for Boost python) 

**General case** : install through your distribution package manager to the system-appropriate directory (/usr, typically)


## Boost

**General case** : install through your distribution package manager to the system-appropriate directory (/usr, typically)

**Special case** : You may want to make the python module as-static-as-possible, that is dynamically linked only to library which
are unaivoidable (like libpython.so or libboost_python.so) but statically linking with all others (other boost libs, rdkit, ...).

See : [Reflections on static or limited runtime dependencies builds](./StaticBuild.md)

## RDKit
**General case** : (no special cases)

Download a RDKit release (https://github.com/rdkit/rdkit/releases).

Unpack somewhere 

Run CMake configuration :

      cd rdkit-Release_2018_09_1
      mkdir RDKit_build
      cd RDKit_build
      
      cmake -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
      -DBoost_USE_STATIC_LIBS=ON \
      -DRDK_BUILD_AVALON_SUPPORT=OFF \
      -DRDK_BUILD_CAIRO_SUPPORT=OFF \
      -DRDK_BUILD_COMPRESSED_SUPPLIERS=OFF \
      -DRDK_BUILD_CONTRIB=OFF \
      -DRDK_BUILD_COORDGEN_SUPPORT=OFF \
      -DRDK_BUILD_CPP_TESTS=OFF \
      -DRDK_BUILD_DESCRIPTORS3D=ON \
      -DRDK_BUILD_FREESASA_SUPPORT=ON \
      -DRDK_BUILD_INCHI_SUPPORT=OFF \
      -DRDK_BUILD_MOLINTERCHANGE_SUPPORT=ON \
      -DRDK_BUILD_PGSQL=OFF \
      -DRDK_BUILD_QT_DEMO=OFF \
      -DRDK_BUILD_QT_SUPPORT=OFF \
      -DRDK_BUILD_RPATH_SUPPORT=OFF \
      -DRDK_BUILD_SLN_SUPPORT=ON \
      -DRDK_BUILD_SWIG_CSHARP_WRAPPER=OFF \
      -DRDK_BUILD_SWIG_JAVA_WRAPPER=OFF \
      -DRDK_BUILD_SWIG_WRAPPERS=OFF \
      -DRDK_BUILD_TEST_GZIP=OFF \
      -DRDK_BUILD_THREADSAFE_SSS=ON \
      -DRDK_COORDGEN_LIBS=MolAlign \
      -DRDK_INSTALL_DEV_COMPONENT=ON \
      -DRDK_INSTALL_DLLS_MSVC=OFF \
      -DRDK_INSTALL_INTREE=ON \
      -DRDK_INSTALL_PYTHON_TESTS=OFF \
      -DRDK_INSTALL_STATIC_LIBS=ON \
      -DRDK_OPTIMIZE_NATIVE=ON \
      -DRDK_PGSQL_STATIC=OFF \
      -DRDK_SWIG_STATIC=OFF \
      -DRDK_TEST_COVERAGE=OFF \
      -DRDK_TEST_MMFF_COMPLIANCE=ON \
      -DRDK_TEST_MULTITHREADED=ON \
      -DRDK_USE_BOOST_REGEX=OFF \
      -DRDK_BUILD_TEST_GZIP=OFF \
      -DRDK_USE_BOOST_SERIALIZATION=OFF \
      -DRDK_USE_FLEXBISON=OFF \
      -DRDK_BUILD_COMPRESSED_SUPPLIERS=OFF \
      -DRDK_USE_STRICT_ROTOR_DEFINITION=ON ..     


Build : 

      make -j4
      make install # dont run as root. Normally does an in tree install but you never know

Some library magic to generate a single static library  :

      cd ../lib
      for libfilename in *.a; do ar -x $libfilename; done
      mkdir obj
      cp *.o obj
      ar cr libRDKitStatic.a *.o
      mkdir staticlib
      cp libRDKitStatic.a staticlib

You may want to change the CMake flags somewhat depending on your needs.
Enabling more flags is generally fine. Disabling is a bit dicey but these are fine (non exhaustive):
- `-DRDK_OPTIMIZE_NATIVE=off`


Releases which have been tested and worked : 
- rdkit-Release_2018_09_1 (no problem detected)
- rdkit_Release_2019_03_1 (need to comment out RDStream `#add_subdirectory(RDStreams)` in rdkit-Release_2019_03_1/Code/CMakeList.txt due to 
static linkage problem) 

## Vc

**General case** : install through your distribution package manager to the system-appropriate directory (/usr, typically)

**Special case** : Build Vc for your target architecture, not -march=native. See Vc documentation.

## libunwind

**General case** : install through your distribution package manager to the system-appropriate directory (/usr, typically)

**Special case** : Build from source with static lib enabled. Not difficult.

## MPI Implementation

**General case** : install through your distribution package manager to the system-appropriate directory (/usr, typically)

OpenMPI has been confirmed to work, but any standard-conforming MPI 3.1 implementation should do the trick.

If needed, here is a ./configure for OpenMPI 4.0.1 that is confirmed to work :

    ./configure --prefix=/home/eliane/local/ --enable-dlopen  --enable-heterogeneous --enable-binaries --enable-shared --enable-static --enable-cxx-exceptions --enable-io-romio

## If needed : External Intel TBB

Generally speaking it is fine to use the shipped version of intel TBB.

However, if you are encountering hangs at start-up; or if you are in an environement where using all CPU
available is undesirable you may : 1/ use the configuration options in SmolDock to configure the number of thread spawned. 2/ 
link your own intel TBB dynamic library instead of the static lib used here.

It may be especially useful if you are running SmolDock concurrently with other intel TBB workload, so as to have a fair work
scheduler. (multiple  executables using the same tbb dynamic libraries "see" the others, and CPU are not oversubscribed)

Build and/or install TBB in the the system-appropriate directory.

# Modifying CMakeList.txt to your local environement

Open CMakeList.txt and go to the part of the file with the following header :

    ###########################################################################
    ######################## BUILD SETUP SECTION ##############################
    ###########################################################################

Then modify the file as described below for each items :

## MPI Implementation

## RDkit

## If needed : External Intel TBB



# Building

    mkdir build
    cd build
    cmake ..
    make -j4 

# Smoke testing

Just run ./smoldock to check if everything works fine : this runs an example docking. You may have to copy data file in the same directory.

The expected output is :
