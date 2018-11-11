cmake_minimum_required(VERSION 3.12)
project(smoldock)

set(CMAKE_CXX_STANDARD 17)



add_library(gsl STATIC IMPORTED)
set_property(TARGET gsl PROPERTY
        IMPORTED_LOCATION /usr/lib/libgsl.a)

add_library(gslcblas STATIC IMPORTED)
set_property(TARGET gslcblas PROPERTY
        IMPORTED_LOCATION /usr/lib/libgslcblas.a)

add_executable(smoldock main.cpp Structures/Structure.cpp Structures/Structure.h Structures/Molecule.cpp Structures/Molecule.h Structures/Protein.cpp Structures/Protein.h Structures/Atom.cpp Structures/Atom.h Structures/Bond.cpp Structures/Bond.h)

target_link_libraries(smoldock gsl gslcblas)