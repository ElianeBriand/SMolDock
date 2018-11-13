cmake_minimum_required(VERSION 3.12)
project(smoldock)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(UnitTests)

add_library(gsl STATIC IMPORTED)
set_property(TARGET gsl PROPERTY
        IMPORTED_LOCATION /usr/lib/libgsl.a)

add_library(gslcblas STATIC IMPORTED)
set_property(TARGET gslcblas PROPERTY
        IMPORTED_LOCATION /usr/lib/libgslcblas.a)

set(SOURCE_LIST Structures/Structure.cpp
        Structures/Structure.h
        Structures/Molecule.cpp
        Structures/Molecule.h
        Structures/Protein.cpp
        Structures/Protein.h
        Structures/Atom.cpp
        Structures/Atom.h
        Structures/Bond.cpp
        Structures/Bond.h
        Utilities/MoleculeTraversal.cpp
        Utilities/MoleculeTraversal.h
        Structures/AminoAcid.cpp
        Structures/AminoAcid.h)

ADD_LIBRARY(smoldock_object OBJECT ${SOURCE_LIST})

add_executable(smoldock main.cpp $<TARGET_OBJECTS:smoldock_object>)

target_link_libraries(smoldock gsl gslcblas)