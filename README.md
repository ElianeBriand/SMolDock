A small protein-ligand docking software.

I'm reimplementing something similar to Autodock Vina because I find the source code to be somewhat
difficult to understand.

Find the documentation on GitHub Pages : https://elianebriand.github.io/SmolDock/

# Build

See the detailed documentation in ./Docs/BuildSmoldock.md or https://elianebriand.github.io/SmolDock/BuildSmoldock.html

In short : Check that you have the required dependencies (see below). Then set the relevant paths in the CMakeLists.txt. Then :

    git submodule update --init --recursive  # init the submodules (tbb)
    mkdir build
    cd build
    cmake ..
    make -j4 # To build all possible frontends and python modules
    # or
    make -j4 smoldock # for the c++ library and test frontend only

Lots of warning but usually nothing blocking.

# Work with it

    cd ./build # if not already in it
    ./smoldock
   
See Frontends/main.cpp for lack of a proper interface yet...

We also have a python interface and various CLI frontend. You want also want to write your own dedicated frontend for your task.

For all that and more, see the documentation in ./Docs or https://elianebriand.github.io/SmolDock/

# Licencing
SmolDock is licenced under GNU GPL version 3 or later.


It includes works from :
- ESBTL (GNU GPLv3)
- RDKit (BSD 3-clauses)
- Autodock Vina (Apache 2.0)
- Ensmallen (BSD 3-clauses)
- Intel TBB (Apache 2.0 - git submodule)

It links with :
- RDKit (BSD 3-clauses)
- Vc (BSD 3-clauses)
- libunwind (MIT licence)
- Eigen (MPL2 licence)
- Boost (Boost software licence)
- Python libraries for the python modules (PSF)
- Any standard-compliant MPI implementation (tested on OpenMPI)

It includes data from :
- GROMACS implementation of Amber99ff (LGPL 2.1 or later)

See COPYING for copyrights and text of these licences.

# Documentation

Checkout the Docs dir. We have code documentation with Doxygen in Docs/html.

Find the documentation online on GitHub Pages : https://elianebriand.github.io/SmolDock/

# Dependencies

## RDKit (linking with pre-built dynamic library)

If you want to be able to export Molecule object to RDKit molecule class, you need to take care that the boost_python 
library used for building RDKit, and SmolDock are exactly the same version (same exact .so version, even minor version number).
Otherwise, the internal conversion will not work. Using a different minor version (.so.1.65.0 vs .so.1.65.1) mostly works,
if you create the appropriate symbolic link, however subtle bug may be introduced (?).



The main author does not rely on pre-built dynamic library shipped by distributions, because of linker
problems that may or may not still exist. See next section for building static libs from source, else :

In the ` RDKIT SETUP ` section of CMakeLists.txt, change :

    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lboost_system -static  ") 

to : 

    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lboost_system ")

And modify the `RDKIT_LINKAGE_LIST` variable, change :
 
    set(RDKIT_LINKAGE_LIST RDKitStatic pthread RDKitStatic RDKitForceField_static RDKitForceFieldHelpers_static)

to whatever is suitable for your installed version of RDKit, probably something like :
 
    set(RDKIT_LINKAGE_LIST RDKit pthread)

or maybe : 

    set(RDKIT_LINKAGE_LIST RDKit pthread RDKitForceField RDKitForceFieldHelpers)


## RDKit (linking with source-built static library)


Download rdkit-Release_2018_09_1 (Garantueed to work, other probably do).

CMake configuration : 

      cd rdkit-Release_2018_09_1
      mkdir RDKit_build
      mkdir rdkit_install
      cd RDKit_build
      
      cmake -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
      -DCMAKE_INSTALL_PREFIX=../rdkit_install
      -DBoost_USE_STATIC_LIBS=OFF \
      -DRDK_BUILD_AVALON_SUPPORT=OFF \
      -DRDK_BUILD_CAIRO_SUPPORT=OFF \
      -DRDK_BUILD_COMPRESSED_SUPPLIERS=OFF \
      -DRDK_BUILD_CONTRIB=ON \
      -DRDK_BUILD_COORDGEN_SUPPORT=ON \
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
      -DRDK_INSTALL_INTREE=OFF \
      -DRDK_INSTALL_PYTHON_TESTS=OFF \
      -DRDK_INSTALL_STATIC_LIBS=OFF \
      -DRDK_OPTIMIZE_NATIVE=ON \
      -DRDK_PGSQL_STATIC=OFF \
      -DRDK_SWIG_STATIC=OFF \
      -DRDK_TEST_COVERAGE=OFF \
      -DRDK_TEST_MMFF_COMPLIANCE=ON \
      -DRDK_TEST_MULTITHREADED=ON \
      -DRDK_USE_BOOST_REGEX=ON \
      -DRDK_USE_BOOST_SERIALIZATION=ON \
      -DRDK_USE_FLEXBISON=OFF \
      -DRDK_USE_STRICT_ROTOR_DEFINITION=ON ..     

Then build : 

      make -j4
      make install

Then library magic :

      cd ../lib
      for libfilename in *.a; do ar -x $libfilename; done
      mkdir obj
      cp *.o obj
      ar cr libRDKitStatic.a *.o
      mkdir staticlib
      cp libRDKitStatic.a staticlib

Edit CMakeLists.txt with the install and source path (search for `RDKIT SETUP`)


For some reason, it is not sufficient to only link against the newly created libRDKit static,
and some additional linkage to the individual .a (like RDKitForceField_static) are necessary as seen in the CMakeLists.txt.
 (Even though theoretically it contains the relevant symbols. A fix/explaination would be welcomed)

If needs be, don't forget you can search where a symbol is defined using

    for filename in *.a; do echo $filename;nm $filename | grep <SYMBOLNAME>; done

The main author is able to create fully statically linked binary using this CMakeLists.txt, in particular
 using `SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lboost_system -static  ")` because she has the relevant libraries
 in static format (libstdc++, libboost_system, libm, libc, libgcc, ...) thanks to a source based distribution. However, if you are
 not in this situation, or have `Undefined reference to ...` linker errors, remove the `-static` from this line : you will still statically
 link to RDKit, but have dependencies on boost, libstdc++ (etc etc) dynamic libraries.

      


## Vc

Edit CMakeLists.txt with your Vc install path (search for `VC SETUP`). You can grab it from GitHub.

## Boost build

If the need arise to have static linking of boost into the python module, this may prove useful to build boost with -fPIC
even for static libraries :


    ./bootstrap.sh --prefix=/home/builder/local/

    ./b2 --prefix=/home/builder/local/ --build-type=complete --build-dir=/home/builder/boost_1_69_0/build/ --layout=tagged cxxflags=" -fPIC " cflags=" -fPIC " --ignore-site-config

## Boost build

If the need arise to have static linking of boost into the python module, this may prove useful to build boost with -fPIC
even for static libraries :


    BOOST_ROOT
    ./bootstrap.sh --prefix=/home/builder/local/

    ./b2 --prefix=/home/builder/local/ --build-type=complete --build-dir=/home/builder/boost_1_69_0/build/ --layout=tagged cxxflags=" -fPIC " cflags=" -fPIC " --ignore-site-config

## Boost build

If the need arise to have static linking of boost into the python module, this may prove useful to build boost with -fPIC
even for static libraries :


    BOOST_ROOT
    ./bootstrap.sh --prefix=/home/builder/local/

    ./b2 --prefix=/home/builder/local/ --build-type=complete --build-dir=/home/builder/boost_1_69_0/build/ --layout=tagged cxxflags=" -fPIC " cflags=" -fPIC " --ignore-site-config