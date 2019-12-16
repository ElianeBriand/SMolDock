A small protein-ligand docking software.

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.com/ElianeBriand/SmolDock.svg?branch=master)](https://travis-ci.com/ElianeBriand/SmolDock)
[![codecov](https://codecov.io/gh/ElianeBriand/SmolDock/branch/master/graph/badge.svg)](https://codecov.io/gh/ElianeBriand/SmolDock)


I'm reimplementing something similar to Autodock Vina because I find the source code to be somewhat
difficult to understand.

Find the documentation on GitHub Pages : https://elianebriand.github.io/SmolDock/

# Build

See the detailed documentation in ./Docs/BuildSmoldock.md or https://elianebriand.github.io/SmolDock/BuildSmoldock.html

Basically, you have to have the following library dependencies installed, including headers, and reachable by the compiler/
build system:

- Intel TBB (Apache 2.0)
- Any standard-compliant MPI implementation (OpenMPI tested to work)
- Vc SIMD (BSD 3-clauses)
- Armadillo (Apache 2.0)
- libunwind (MIT licence)
- Eigen (MPL2 licence)
- Boost (Boost software licence)
- Python libraries for the python modules (PSF)


These can be installed through your package manager (you might need the \*-devel package). For example, for Ubuntu bionic you would need something like :

```
apt-get install g++-8 g++-8-multilib libhwloc-plugins openmpi-bin libopenmpi-dev libboost-all-dev lcov zlib1g-dev libzip-dev libbz2-dev lzma liblzma-dev libtbb-dev libopenmpi-dev vc-dev libzstd-dev xz-utils libunwind-dev libarmadillo-dev libeigen3-dev libpython3.7-dev python3.7-dev wget libsqlite3-dev
```

(See alos .travis.yml, as the CI script generally is able to build and run the software)

You also need to have RDKit build with a specific set of parameters. It is suggested to download the latest release and 
adapt the following instruction :

      tar xf ReleaseXXX
      cd rdkit-Release_XXX
      mkdir rdkit_build
      mkdir rdkit_install
      cd rdkit_build
      
      cmake -DCMAKE_INSTALL_PREFIX=../rdkit_install \
      -DRDK_INSTALL_INTREE=OFF \
      -DBUILD_TESTING=OFF \
      -DCMAKE_BUILD_TYPE=Release \
      -DRDK_BUILD_AVALON_SUPPORT=ON \
      -DRDK_BUILD_CAIRO_SUPPORT=OFF \
      -DRDK_BUILD_COMPRESSED_SUPPLIERS=ON \
      -DRDK_BUILD_CONTRIB=ON \
      -DRDK_BUILD_COORDGEN_SUPPORT=ON \
      -DRDK_BUILD_CPP_TESTS=OFF \
      -DRDK_BUILD_DESCRIPTORS3D=ON \
      -DRDK_BUILD_FREESASA_SUPPORT=ON \
      -DRDK_BUILD_INCHI_SUPPORT=OFF \
      -DRDK_BUILD_MINIMAL_LIB=OFF \
      -DRDK_BUILD_MOLINTERCHANGE_SUPPORT=ON \
      -DRDK_BUILD_PGSQL=OFF \
      -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
      -DRDK_BUILD_QT_DEMO=OFF \
      -DRDK_BUILD_QT_SUPPORT=OFF \
      -DRDK_BUILD_RPATH_SUPPORT=OFF \
      -DRDK_BUILD_SLN_SUPPORT=ON \
      -DRDK_BUILD_SWIG_CSHARP_WRAPPER=OFF    \
      -DRDK_BUILD_SWIG_JAVA_WRAPPER=OFF      \
      -DRDK_BUILD_SWIG_WRAPPERS=OFF          \
      -DRDK_BUILD_TEST_GZIP=OFF              \
      -DRDK_BUILD_THREADSAFE_SSS=ON          \
      -DRDK_BUILD_YAEHMOP_SUPPORT=OFF         \
      -DRDK_INSTALL_DEV_COMPONENT=ON         \
      -DRDK_INSTALL_DLLS_MSVC=OFF            \
      -DRDK_INSTALL_PYTHON_TESTS=OFF         \
      -DRDK_INSTALL_STATIC_LIBS=OFF          \
      -DRDK_OPTIMIZE_NATIVE=ON               \
      -DRDK_PGSQL_STATIC=OFF                 \
      -DRDK_SWIG_STATIC=OFF                  \
      -DRDK_TEST_COVERAGE=OFF                \
      -DRDK_TEST_MMFF_COMPLIANCE=OFF         \
      -DRDK_TEST_MULTITHREADED=OFF           \
      -DRDK_USE_BOOST_IOSTREAMS=ON           \
      -DRDK_USE_BOOST_REGEX=ON               \
      -DRDK_USE_BOOST_SERIALIZATION=ON       \
      -DRDK_USE_FLEXBISON=OFF                \
      -DRDK_USE_STRICT_ROTOR_DEFINITION=ON    \
      -DRDK_USE_URF=ON                       \
      ..
      
      make -j4
      # Check if that previous command completed successfully
      # If not, try running it again
      # If still not, you might need to install dependencies for RDKit
      # You're on your own here, look at the error message
      make install
 

Then you will need to edit the RDKIT_ROOT path in the CMakeList.txt. It is at the top of the file, in the build configuration
section. You'll probably need to change the other variable in that section : follow the instruction there.

Then, clone SmolDock and init the submodules:

    git clone https://github.com/ElianeBriand/SmolDock.git
    git submodule update --init --recursive

Then build:

    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
    make -j4


# Use

    cd ./build # if not already in it
    ./smoldock
   
At this moment, there is no real front-end, so modify Frontends/main.cpp then recompile.

Also: see the documentation in ./Docs or https://elianebriand.github.io/SmolDock/

# Licencing

SmolDock is licenced under GNU GPL version 3 or later.

It includes works from :
- ESBTL (GNU GPLv3)
- RDKit (BSD 3-clauses)
- Autodock Vina (Apache 2.0)
- Ensmallen (BSD 3-clauses)


It links with :
- RDKit (BSD 3-clauses)
- Vc (BSD 3-clauses)
- libunwind (MIT licence)
- Eigen (MPL2 licence)
- Boost (Boost software licence)
- Python libraries for the python modules (PSF)
- Intel TBB (Apache 2.0)
- Armadillo (Apache 2.0)
- Any standard-compliant MPI implementation (tested on OpenMPI)
- PDBPC (GPLv3)

It includes data from :
- GROMACS implementation of Amber99ff (LGPL 2.1 or later)

See COPYING for copyrights and text of these licences.

# Documentation

Checkout the Docs dir. We have code documentation with Doxygen in Docs/html.

Find the documentation online on GitHub Pages : https://elianebriand.github.io/SmolDock/

