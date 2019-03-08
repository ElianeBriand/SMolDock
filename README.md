A small protein-ligand docking software.

I'm reimplementing something similar to Autodock Vina because I find the source code to be somewhat
difficult to understand.

Find the documentation on GitHub Pages : https://elianebriand.github.io/SmolDock/

# Build

Check that you have the required dependencies (see below), then : 

    cmake .
    make -j4

# Work with it

    cd ./cmake-build-XXX
    ./smoldock
   
See Frontends/main.cpp for lack of a proper interface yet...


# Licencing
SmolDock is licenced under GNU GPL version 3 or later.


It includes works from :

- ESBTL (GNU GPLv3)
- Autodock Vina (Apache 2.0)
- RDKit (BSD 3-clauses)
- Ensmallen (BSD 3-clauses)
- Eigen (MPL2 licence)

It links with
- RDKit (BSD 3-clauses)
- Vc (BSD 3-clauses)

It includes data from :
- GROMACS implementation of Amber99ff (LGPL 2.1 or later)

See COPYING for copyrights and text of these licences.

# Documentation

Checkout the Docs dir. We have code documentation with Doxygen in Docs/html, and sphinx-based more general HOWTO
in Docs/build/html.

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
      cd RDKit_build
      
      cmake -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
      -DBoost_USE_STATIC_LIBS=ON \
      -DRDK_BUILD_AVALON_SUPPORT=OFF \
      -DRDK_BUILD_CAIRO_SUPPORT=OFF \
      -DRDK_BUILD_COMPRESSED_SUPPLIERS=OFF \
      -DRDK_BUILD_CONTRIB=OFF \
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
      -DRDK_USE_BOOST_SERIALIZATION=OFF \
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


    BOOST_ROOT
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