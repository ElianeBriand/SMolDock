cmake_minimum_required(VERSION 3.12)

################ BUILD CONFIGURATION SECTION     ############################

# Point this to the rdkit install directory (See Readme)
set(RDKIT_ROOT /home/eliane/Projects/rdkit_releases/rdkit-Release_2019_09_1/rdkit_install/)


# Command to run : mpic++ -showme:compile. Only add the argument part of the command
set(MPI_COMPILE_FLAGS -pthread -Wl,-rpath -Wl,/usr/lib64 -Wl,--enable-new-dtags -L/usr/lib64 -lmpi_cxx -lmpi)

#Command to run  : mpic++  -showme:link
set(MPI_LINK_FLAGS  -pthread -Wl,-rpath -Wl,/usr/lib64 -Wl,--enable-new-dtags -L/usr/lib64 -lmpi_cxx -lmpi)

# Version of python installed with header and libraries
set(PY_VERSION 3.6)

################ END BUILD CONFIGURATION SECTION ############################


set(CMAKE_CXX_STANDARD 17)

#set(BOOST_ROOT /home/briand/local)
#set(PYTHON_INCLUDE_DIR /home/briand/local/include/python3.6m)
#set(PYTHON_LIBRARY /home/briand/local/lib)
#SET(CMAKE_CXX_COMPILER /home/briand/local/bin/g++)

# where is the target environment
#SET(CMAKE_FIND_ROOT_PATH  /home/briand/local/)
#SET(CMAKE_PREFIX_PATH /home/briand/local/)

SET(CMAKE_CXX_FLAGS "-Wall")
SET(CMAKE_C_COMPILER  "/usr/bin/g++")
SET(CMAKE_CXX_COMPILER "/usr/bin/g++")
SET(CMAKE_CXX_FLAGS_DEBUG "-g ")
SET(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -march=native -g  ")

set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
            STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
            "Debug" "Release" "RelWithDebInfo")
endif()


project(smoldock)





###########################################################################
######################## BUILD SETUP SECTION ##############################
###########################################################################



## MPI #######################
# Edit this if you have a specific MPI setup
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

## END MPI ###################

## Boost #######################
# Edit this if you have a specific Boost setup
find_package(Boost REQUIRED regex)
include_directories(${Boost_INCLUDE_DIRS})
#add_compile_definitions(BOOST_ALLOW_DEPRECATED_HEADERS)
## END Boost ###################

## Python #######################
# Edit version and other setup if needed

find_package(PythonLibs ${PY_VERSION} REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
## END Python ###################

######## VC SETUP ######
# Adjust this according to where your Vc build is :
#set(VC_ROOT /usr/)
#
#include_directories(${VC_ROOT}/include/)
#LINK_DIRECTORIES(${VC_ROOT}/lib/)
######## END OF VC  ######

######## EIGEN SETUP ######
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
#set(EIGEN_ROOT /usr/)
#include_directories(${EIGEN_ROOT}/include)
######## END OF EIGEN  ######

######## RDKIT SETUP ######
# Adjust this according to where your RDKit build is :
#set(RDKIT_ROOT /home/eliane/Projects/RDKit_Git/rdkit-Release_2018_09_1)

# This directory should contain ./Code ./lib and ./lib/staticlib as setup in the build guide
include(${RDKIT_ROOT}/lib/cmake/rdkit/rdkit-config.cmake)
include_directories(${RDKIT_ROOT}/include/)
include_directories(${RDKIT_ROOT}/include/rdkit)
LINK_DIRECTORIES(${RDKIT_ROOT}/lib)
######## END OF RDKIT ######



######## TBB SETUP ######
find_package( Threads )
#add_subdirectory(Dependencies/tbb)
#include_directories(Dependencies/tbb/include)
######## END OF TBB ######



# Comment out the following for quieter runs (generally speaking, do not affect runtime either way)
add_compile_definitions(SMOLDOCK_VERBOSE_DEBUG)

# For Intel VTune profiler tools, set this to true (necessary if building a static executable that can be profiled)
set(SMOLDOCK_INTEL_VTUNE_EXPORT_SYMBOL FALSE)

###########################################################################
######################## END  BUILD SETUP SECTION #########################
###########################################################################


include_directories(Dependencies)
include_directories(.)



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
        Structures/AminoAcid.cpp
        Structures/AminoAcid.h
        Engines/AbstractDockingEngine.h
        Engines/MDStyleDockingEngine.cpp
        Engines/MDStyleDockingEngine.h
        Engines/VinaCompatibleDockingEngine.cpp
        Engines/VinaCompatibleDockingEngine.h
        Structures/Results/DockingResult.cpp
        Structures/Results/DockingResult.h
        Utilities/DockingResultPrinter.cpp
        Utilities/DockingResultPrinter.h
        Engines/ConformerDockingEngine.cpp
        Engines/ConformerDockingEngine.h
        Engines/Internals/iConformer.h
        Utilities/TimingsLog.h
        Engines/ScoringFunctions/VinaLikeRigid.cpp
        Engines/ScoringFunctions/VinaLikeRigid.h
        Engines/Internals/iProtein.h
        Engines/Internals/iTransform.h
        Utilities/PDBWriter.cpp
        Utilities/PDBWriter.h Utilities/PDBLigandUtils.h
        Utilities/PDBLigandUtils.cpp
        Engines/LocalOptimizers/GradientDescentLineSearch.cpp
        Engines/LocalOptimizers/GradientDescentLineSearch.h
        Utilities/IntermediateConformerCollector.cpp
        Utilities/IntermediateConformerCollector.h
        Utilities/ReScorer.cpp
        Utilities/ReScorer.h
        Structures/Common/ResiduePropertiesAssignation.cpp
        Structures/Common/ResiduePropertiesAssignation.h
        Structures/Common/PDBResiduePropertiesTable.cpp
        Structures/Common/PDBResiduePropertiesTable.h
        Structures/InputModifiers/VinaCompatibility.cpp
        Structures/InputModifiers/VinaCompatibility.h
        Structures/InputModifiers/InputModifierInterface.h
        Engines/ScoringFunctions/ScoringFunctionInterface.h
        Engines/LocalOptimizers/L_BFGS.cpp
        Engines/LocalOptimizers/L_BFGS.h
        Engines/LocalOptimizers/OptimizerInterface.h
        Engines/Internals/InternalsUtilityFunctions.h
        Engines/LocalOptimizers/LocalOptimizerUtilityFunctions.h
        Engines/GlobalHeuristics/RandomRestart.cpp
        Engines/GlobalHeuristics/RandomRestart.h
        Engines/GlobalHeuristics/HeuristicInterface.h
        Engines/ScoringFunctions/ScoringFunctionFactory.cpp
        Engines/ScoringFunctions/ScoringFunctionFactory.h
        Engines/LocalOptimizers/OptimizerFactory.cpp
        Engines/LocalOptimizers/OptimizerFactory.h
        Engines/GlobalHeuristics/HeuristicFactory.cpp
        Engines/GlobalHeuristics/HeuristicFactory.h
        Engines/PoseRefiner.cpp
        Engines/PoseRefiner.h
        Engines/GlobalHeuristics/OnlyLocal.cpp
        Engines/GlobalHeuristics/OnlyLocal.h
        Engines/GlobalHeuristics/IteratedLocalSearch.cpp
        Engines/GlobalHeuristics/IteratedLocalSearch.h
        Engines/GlobalHeuristics/SimulatedAnnealing.cpp
        Engines/GlobalHeuristics/SimulatedAnnealing.h
        Engines/GlobalHeuristics/Utilities/MetropolisHastings.cpp
        Engines/GlobalHeuristics/Utilities/MetropolisHastings.h
        Utilities/Version.cpp
        Utilities/Version.h
        Engines/ScoringFunctions/VinaLike.cpp
        Engines/ScoringFunctions/VinaLike.h
        Engines/ScoringFunctions/VinaLikeCommon.h
        Engines/GlobalHeuristics/DifferentialEvolution.cpp
        Engines/GlobalHeuristics/DifferentialEvolution.h
        Engines/GlobalHeuristics/Evolution.cpp
        Engines/GlobalHeuristics/Evolution.h
        Engines/ScoringFunctions/VinaLikeCovalentReversible.cpp
        Engines/ScoringFunctions/VinaLikeCovalentReversible.h
        Structures/InputModifiers/RotatableBondRemover.cpp
        Structures/InputModifiers/RotatableBondRemover.h
        Engines/ScoringFunctions/VinaLikeCovalentReversible.cpp
        Engines/ScoringFunctions/VinaLikeCovalentReversible.h
        Utilities/Calibration/Calibrator.cpp
        Utilities/Calibration/Calibrator.h
        Utilities/CSVReader.cpp Utilities/CSVReader.h
        Utilities/SMARTSMatcher.cpp
        Utilities/SMARTSMatcher.h
        Utilities/AdvancedErrorHandling.cpp
        Utilities/AdvancedErrorHandling.h)

set(MPI_SOURCE_LIST Utilities/Calibration/MPICalibratorDirector.cpp
        Utilities/Calibration/MPICalibratorDirector.h
        Utilities/Calibration/MPICalibratorNode.cpp
        Utilities/Calibration/MPICalibratorNode.h
        Utilities/Calibration/MPICalibratorCommon.h
        )

set(PY_EXPORT_SOURCE_LIST Frontends/Python/PyStructures.cpp
        Frontends/Python/PyStructures.h
        Frontends/Python/PySTLWrapper.cpp
        Frontends/Python/PySTLWrapper.h
        Frontends/Python/PyUtilities.cpp
        Frontends/Python/PyUtilities.h
        Frontends/Python/PyEngine.cpp
        Frontends/Python/PyEngine.h Frontends/FrontendCommon.h)




set(BPRINTER_SOURCE_LIST
        Dependencies/bprinter/table_printer.h
        Dependencies/bprinter/impl/table_printer.tpp.h
        Dependencies/bprinter/table_printer.cpp
        )

add_compile_definitions(BOOST_ENABLE_ASSERT_HANDLER)


ADD_LIBRARY(libsmoldock SHARED ${SOURCE_LIST} ${BPRINTER_SOURCE_LIST})
target_compile_definitions(libsmoldock PUBLIC ARMA_DONT_USE_WRAPPER)
add_dependencies(libsmoldock rdkit_base)
target_compile_definitions(libsmoldock PUBLIC BOOST_LOG_DYN_LINK)


add_executable(smoldock Frontends/main.cpp)
target_link_libraries(smoldock libsmoldock)
target_link_libraries(smoldock boost_log_setup boost_log boost_thread)
target_link_libraries(smoldock tbb)
target_link_libraries(smoldock rdkit_base avalon_clib AvalonLib freesasa_clib FreeSASALib maeparser )
target_link_libraries(smoldock coordgen RingDecomposerLib URFLib RDGeneral RDStreams DataStructs )
target_link_libraries(smoldock RDGeometryLib Alignment EigenSolvers Optimizer ForceField DistGeometry)
target_link_libraries(smoldock Catalogs GraphMol Depictor SmilesParse FileParsers SubstructMatch ChemReactions ChemTransforms Subgraphs)
target_link_libraries(smoldock FilterCatalog FragCatalog Descriptors Fingerprints PartialCharges MolTransforms ForceFieldHelpers)
target_link_libraries(smoldock DistGeomHelpers MolAlign MolChemicalFeatures ShapeHelpers MolCatalog MolDraw2D FMCS MolHash MMPA)
target_link_libraries(smoldock StructChecker ReducedGraphs Trajectory SubstructLibrary RGroupDecomposition MolInterchange SLNParse)
target_link_libraries(smoldock MolStandardize SimDivPickers hc InfoTheory ChemicalFeatures ConformerParser)
target_link_libraries(smoldock unwind)
target_link_libraries(smoldock Vc)
target_link_libraries(smoldock lapack blas)
target_link_libraries(smoldock pthread)


if(SMOLDOCK_INTEL_VTUNE_EXPORT_SYMBOL)
    target_link_options(smoldock PRIVATE "-u_init" "-usetenv" "-ugetenv" "-u__errno_location" "-upthread_create" "-upthread_getattr_np")
    target_link_options(smoldock PRIVATE "-upthread_self" "-upthread_getspecific" "-upthread_setspecific")
    target_link_options(smoldock PRIVATE "-upthread_getattr_np" "-upthread_attr_destroy" "-upthread_attr_setstack" "-upthread_attr_getstack")
    target_link_options(smoldock PRIVATE "-upthread_attr_getstacksize" "-upthread_attr_setstacksize" "-upthread_cancel" "-u_pthread_cleanup_push")
    target_link_options(smoldock PRIVATE "-u_pthread_cleanup_pop" "-upthread_attr_setstacksize" "-upthread_cancel" "-u_pthread_cleanup_push")
    target_link_options(smoldock PRIVATE "-upthread_mutex_lock" "-upthread_mutex_trylock" "-upthread_spin_lock" "-upthread_spin_trylock")
endif()


add_executable(goldenRecordGenerator Frontends/goldenRecordGenerator.cpp)
target_link_libraries(goldenRecordGenerator libsmoldock)
target_link_libraries(goldenRecordGenerator boost_log_setup boost_log boost_thread)
target_link_libraries(goldenRecordGenerator tbb)
target_link_libraries(goldenRecordGenerator rdkit_base avalon_clib AvalonLib freesasa_clib FreeSASALib maeparser )
target_link_libraries(goldenRecordGenerator coordgen RingDecomposerLib URFLib RDGeneral RDStreams DataStructs )
target_link_libraries(goldenRecordGenerator RDGeometryLib Alignment EigenSolvers Optimizer ForceField DistGeometry)
target_link_libraries(goldenRecordGenerator Catalogs GraphMol Depictor SmilesParse FileParsers SubstructMatch ChemReactions ChemTransforms Subgraphs)
target_link_libraries(goldenRecordGenerator FilterCatalog FragCatalog Descriptors Fingerprints PartialCharges MolTransforms ForceFieldHelpers)
target_link_libraries(goldenRecordGenerator DistGeomHelpers MolAlign MolChemicalFeatures ShapeHelpers MolCatalog MolDraw2D FMCS MolHash MMPA)
target_link_libraries(goldenRecordGenerator StructChecker ReducedGraphs Trajectory SubstructLibrary RGroupDecomposition MolInterchange SLNParse)
target_link_libraries(goldenRecordGenerator MolStandardize SimDivPickers hc InfoTheory ChemicalFeatures ConformerParser)
target_link_libraries(goldenRecordGenerator unwind)
target_link_libraries(goldenRecordGenerator Vc)
target_link_libraries(goldenRecordGenerator lapack blas)
target_link_libraries(goldenRecordGenerator pthread)


add_executable(mpicalibrator Frontends/mpicalibrator.cpp ${MPI_SOURCE_LIST})
set_target_properties(mpicalibrator PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS} ")
set_target_properties(mpicalibrator PROPERTIES LINK_FLAGS ${MPI_LINK_FLAGS})
target_link_libraries(mpicalibrator libsmoldock)
target_link_libraries(mpicalibrator boost_log_setup boost_log boost_thread)
target_link_libraries(mpicalibrator tbb)
target_link_libraries(mpicalibrator rdkit_base avalon_clib AvalonLib freesasa_clib FreeSASALib maeparser )
target_link_libraries(mpicalibrator coordgen RingDecomposerLib URFLib RDGeneral RDStreams DataStructs )
target_link_libraries(mpicalibrator RDGeometryLib Alignment EigenSolvers Optimizer ForceField DistGeometry)
target_link_libraries(mpicalibrator Catalogs GraphMol Depictor SmilesParse FileParsers SubstructMatch ChemReactions ChemTransforms Subgraphs)
target_link_libraries(mpicalibrator FilterCatalog FragCatalog Descriptors Fingerprints PartialCharges MolTransforms ForceFieldHelpers)
target_link_libraries(mpicalibrator DistGeomHelpers MolAlign MolChemicalFeatures ShapeHelpers MolCatalog MolDraw2D FMCS MolHash MMPA)
target_link_libraries(mpicalibrator StructChecker ReducedGraphs Trajectory SubstructLibrary RGroupDecomposition MolInterchange SLNParse)
target_link_libraries(mpicalibrator MolStandardize SimDivPickers hc InfoTheory ChemicalFeatures ConformerParser)
target_link_libraries(mpicalibrator unwind)
target_link_libraries(mpicalibrator Vc)
target_link_libraries(mpicalibrator lapack blas)
target_link_libraries(mpicalibrator boost_mpi mpi_cxx mpi)
target_link_libraries(mpicalibrator pthread)




add_library(PySmolDock MODULE Frontends/SmoldockPythonBindings.cpp ${PY_EXPORT_SOURCE_LIST})
set_target_properties(PySmolDock PROPERTIES COMPILE_FLAGS " -fPIC ")
set_target_properties(PySmolDock PROPERTIES LINK_FLAGS " -shared ")
set_target_properties(PySmolDock PROPERTIES PREFIX "")
set_target_properties(PySmolDock PROPERTIES OUTPUT_NAME "PySmolDock")
target_compile_definitions(PySmolDock PUBLIC BOOST_LOG_DYN_LINK)
target_compile_definitions(PySmolDock PUBLIC ARMA_DONT_USE_WRAPPER)
target_link_libraries(PySmolDock libsmoldock)
target_link_libraries(PySmolDock boost_log_setup boost_log boost_thread)
target_link_libraries(PySmolDock tbb)
target_link_libraries(PySmolDock rdkit_base avalon_clib AvalonLib freesasa_clib FreeSASALib maeparser )
target_link_libraries(PySmolDock coordgen RingDecomposerLib URFLib RDGeneral RDStreams DataStructs )
target_link_libraries(PySmolDock RDGeometryLib Alignment EigenSolvers Optimizer ForceField DistGeometry)
target_link_libraries(PySmolDock Catalogs GraphMol Depictor SmilesParse FileParsers SubstructMatch ChemReactions ChemTransforms Subgraphs)
target_link_libraries(PySmolDock FilterCatalog FragCatalog Descriptors Fingerprints PartialCharges MolTransforms ForceFieldHelpers)
target_link_libraries(PySmolDock DistGeomHelpers MolAlign MolChemicalFeatures ShapeHelpers MolCatalog MolDraw2D FMCS MolHash MMPA)
target_link_libraries(PySmolDock StructChecker ReducedGraphs Trajectory SubstructLibrary RGroupDecomposition MolInterchange SLNParse)
target_link_libraries(PySmolDock MolStandardize SimDivPickers hc InfoTheory ChemicalFeatures ConformerParser)
target_link_libraries(PySmolDock unwind)
target_link_libraries(PySmolDock Vc)
target_link_libraries(PySmolDock lapack blas)
target_link_libraries(PySmolDock pthread)

IF(CMAKE_BUILD_TYPE MATCHES Debug)
    target_compile_definitions(libsmoldock PUBLIC SMOLDOCK_VERBOSE_DEBUG)
    target_compile_definitions(smoldock PUBLIC SMOLDOCK_VERBOSE_DEBUG)
    target_compile_definitions(PySmolDock PUBLIC SMOLDOCK_VERBOSE_DEBUG)
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)




###############################################################################################################
############                                  UNIT TESTING
###############################################################################################################
enable_testing()

set(TEST_SOURCE_LIST
        UnitTests/main_test_suite.cpp
        UnitTests/engine_test_suite.cpp)


add_executable(complete_test_suite ${TEST_SOURCE_LIST})

set_target_properties(complete_test_suite PROPERTIES LINK_FLAGS "   ")
target_link_libraries(complete_test_suite libsmoldock)
target_link_libraries(complete_test_suite boost_log_setup boost_log boost_thread)
target_link_libraries(complete_test_suite tbb)
target_link_libraries(complete_test_suite rdkit_base avalon_clib AvalonLib freesasa_clib FreeSASALib maeparser )
target_link_libraries(complete_test_suite coordgen RingDecomposerLib URFLib RDGeneral RDStreams DataStructs )
target_link_libraries(complete_test_suite RDGeometryLib Alignment EigenSolvers Optimizer ForceField DistGeometry)
target_link_libraries(complete_test_suite Catalogs GraphMol Depictor SmilesParse FileParsers SubstructMatch ChemReactions ChemTransforms Subgraphs)
target_link_libraries(complete_test_suite FilterCatalog FragCatalog Descriptors Fingerprints PartialCharges MolTransforms ForceFieldHelpers)
target_link_libraries(complete_test_suite DistGeomHelpers MolAlign MolChemicalFeatures ShapeHelpers MolCatalog MolDraw2D FMCS MolHash MMPA)
target_link_libraries(complete_test_suite StructChecker ReducedGraphs Trajectory SubstructLibrary RGroupDecomposition MolInterchange SLNParse)
target_link_libraries(complete_test_suite MolStandardize SimDivPickers hc InfoTheory ChemicalFeatures ConformerParser)
target_link_libraries(complete_test_suite unwind)
target_link_libraries(complete_test_suite Vc)
target_link_libraries(complete_test_suite lapack blas)
target_link_libraries(complete_test_suite pthread)
target_link_libraries(complete_test_suite ${Boost_LIBRARIES} boost_unit_test_framework )
#target_link_libraries(complete_test_suite Vc)


