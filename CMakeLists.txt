cmake_minimum_required(VERSION 3.12)
project(smoldock)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})

if (COVERAGE)
    SET(GCC_COVERAGE_COMPILE_FLAGS "-fprofile-arcs -ftest-coverage")
endif ()

if (SANITIZER)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer ")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
endif ()

######## RDKIT ######
# Adjust this according to where your RDKit build is :
set(RDKIT_ROOT /home/eliane/Projects/RDKit/)

include_directories(/home/eliane/Projects/rdkit_src/rdkit/Code/)

LINK_DIRECTORIES(${RDKIT_ROOT}/lib)
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lboost_system ")
######## END OF RDKIT ######


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


target_link_libraries(smoldock ${Boost_LIBRARIES})

# RDKit lib*.so
TARGET_LINK_LIBRARIES(smoldock RDKitAlignment RDKitForceFieldHelpers RDKitRDBoost RDKitCatalogs RDKitForceField RDKitRDGeneral RDKitChemicalFeatures RDKitFragCatalog RDKitRDGeometryLib RDKitChemReactions RDKitGraphMol RDKitReducedGraphs RDKitChemTransforms RDKithc RDKitRGroupDecomposition RDKitDataStructs RDKitInfoTheory RDKitShapeHelpers RDKitDepictor RDKitMMPA RDKitSimDivPickers RDKitDescriptors RDKitMolAlign RDKitSLNParse RDKitDistGeometry RDKitMolCatalog RDKitSmilesParse RDKitDistGeomHelpers RDKitMolChemicalFeatures RDKitStructChecker RDKitEigenSolvers RDKitMolDraw2D RDKitSubgraphs RDKitFileParsers RDKitMolHash RDKitSubstructLibrary RDKitFilterCatalog RDKitMolTransforms RDKitSubstructMatch RDKitFingerprints RDKitOptimizer RDKitTrajectory RDKitFMCS RDKitPartialCharges)


target_link_libraries(smoldock gsl gslcblas)