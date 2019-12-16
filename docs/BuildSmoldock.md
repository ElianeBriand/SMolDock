
The build environement for SmolDock is  somewhat cumbersome, as there are a relative lot of dependencies, but apart from
RDkit, the installation of these should be relatively straightforward

# Dependencies likely to be found in your package manager

You need the following library installed, with headers too. You may need to install the *-devel package for each of them
depending on your linux distribution. On Windows, you may have to collect/build them manually (which may be a lot of work).

- Intel TBB (Apache 2.0)
- Any standard-compliant MPI implementation (OpenMPI tested to work)
- Vc SIMD (BSD 3-clauses)
- libunwind (MIT licence)
- Eigen (MPL2 licence)
- Boost (Boost software licence)
- Python libraries for the python modules (PSF)

# Building RDKit

Probably the most difficult part. You also need to have RDKit build with a specific set of parameters. It is suggested
to download the latest release and  adapt the following instruction: (see below in case of problem)
 
 
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

## Possible problems while building RDKit

If the error message mentions missing library or symbols, you may have to install the devel version of the mentionned libraries.

# Editing CMakeList.txt

Look at the top of the CMakeList.txt, and change the variable to suit your system.
