FROM ubuntu:18.04

LABEL maintainer="Eliane Briand <eliane@br.iand.fr>"

# disable apt-get questions
ARG DEBIAN_FRONTEND=noninteractive

ARG cmake_latest_deps="apt-transport-https ca-certificates gnupg software-properties-common wget sudo"
ARG utilities_deps="cmake sudo build-essential gfortran g++ ca-certificates"
ARG build_deps="libboost-all-dev libopenmpi-dev libhwloc-plugins lcov zlib1g-dev libzip-dev libbz2-dev lzma liblzma-dev libtbb-dev vc-dev libzstd-dev xz-utils libunwind-dev libarmadillo-dev libeigen3-dev libpython3.7-dev python3.7-dev libsqlite3-dev"

# install libraries 
RUN apt-get -yq update \
 && apt-get -yq install --no-install-recommends $cmake_latest_deps \
 && wget -qO - https://apt.kitware.com/keys/kitware-archive-latest.asc | sudo apt-key add - \
 && sudo apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main' \
 && sudo apt-get update \
 && apt-get -yq install --no-install-recommends $utilities_deps $build_deps

RUN mkdir -p /tmp/smoldock_build
ADD . /tmp/smoldock_build
RUN export PYTHON_EXECUTABLE=`which python3` && export PYTHON_INCLUDE_DIR=$(python3 -c "from distutils.sysconfig import get_python_inc; print(get_python_inc())") && export PYTHON_LIBRARY=$(python3 -c "import distutils.sysconfig as sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
RUN wget https://github.com/VcDevel/Vc/releases/download/1.4.1/Vc-1.4.1.tar.gz \
    && tar xf Vc-1.4.1.tar.gz \
    && cd Vc-1.4.1 \
    && mkdir Vc_build \
    && cd Vc_build \
    && cmake -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=/usr/ -DUSE_LIBMVEC=OFF .. \
    && make -j "$(nproc)" \
    && sudo make install \
    && cd ../..
RUN wget https://github.com/rdkit/rdkit/archive/Release_2019_09_2.tar.gz \
    && mkdir rdkit_install \
    && tar xf Release_2019_09_2.tar.gz \
    && cd rdkit-Release_2019_09_2 \
    && mkdir rdkit_build \
    && cd rdkit_build \
    && cmake -DCMAKE_INSTALL_PREFIX=../../rdkit_install \
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
    -DRDK_USE_URF=ON  \
    .. \
    && make -j "$(nproc)" \
    && make install \
    && cd ../..
RUN export LD_LIBRARY_PATH=`pwd`/rdkit_install/lib/:$LD_LIBRARY_PATH \
    && mkdir build \
    && cd build \
    && cmake -DRDKIT_ROOT=`pwd`/../rdkit_install /tmp/smoldock_build \
    && make -j "$(nproc)" \
    && cd .. \
    && sleep 10


