cmake_minimum_required(VERSION 3.12)

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

# Use mpic++ -showme (or your MPI implementation variant) to get the corresponding compile and link flags

# Command to run : mpic++ -showme:compile
set(MPI_COMPILE_FLAGS -pthread)

#Command to run  : mpic++ -static -showme:link
# you may want to add these afterward if using openMPI : "mpi open-rte open-pal m dl z rt util"
set(MPI_LINK_FLAGS  -pthread -Wl,-rpath -Wl,/usr/lib64 -Wl,--enable-new-dtags -L/usr/lib64 -lmpi_cxx -lmpi)
## END MPI ###################

## Boost #######################
# Edit this if you have a specific Boost setup
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
#add_compile_definitions(BOOST_ALLOW_DEPRECATED_HEADERS)
## END Boost ###################

## Python #######################
# Edit version and other setup if needed
set(PY_VERSION 3.7)
find_package(PythonLibs ${PY_VERSION} REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
## END Python ###################

######## VC SETUP ######
# Adjust this according to where your Vc build is :
set(VC_ROOT /usr/)

include_directories(${VC_ROOT}/include/)
LINK_DIRECTORIES(${VC_ROOT}/lib/)
######## END OF VC  ######

######## EIGEN SETUP ######
# Adjust this according to where your Vc build is :
set(EIGEN_ROOT /usr/)

include_directories(${EIGEN_ROOT}/include)
######## END OF EIGEN  ######

######## RDKIT SETUP ######
# Adjust this according to where your RDKit build is :
#set(RDKIT_ROOT /home/eliane/Projects/RDKit_Git/rdkit-Release_2018_09_1)
set(RDKIT_ROOT /home/eliane/Projects/rdkit_releases/rdkit-Release_2019_03_1/)

# This directory should contain ./Code ./lib and ./lib/staticlib as setup in the build guide

include_directories(${RDKIT_ROOT}/Code/)
LINK_DIRECTORIES(${RDKIT_ROOT}/lib/staticlib)
LINK_DIRECTORIES(${RDKIT_ROOT}/lib)
set(RDKIT_LINKAGE_LIST RDKitStatic RDKitStatic)
#set(RDKIT_LINKAGE_LIST RDKitStatic RDKitStatic RDKitForceFieldHelpers_static RDKitForceField_static)
######## END OF RDKIT ######


######## TBB SETUP ######
add_subdirectory(Dependencies/tbb)
include_directories(Dependencies/tbb/include)
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




add_library(bprinter STATIC
        Dependencies/bprinter/table_printer.h
        Dependencies/bprinter/impl/table_printer.tpp.h
        Dependencies/bprinter/table_printer.cpp
        )

add_compile_definitions(BOOST_ENABLE_ASSERT_HANDLER)


ADD_LIBRARY(smoldock_object OBJECT ${SOURCE_LIST})
set_target_properties(smoldock_object PROPERTIES LINK_FLAGS " -static  ")
target_compile_definitions(smoldock_object PUBLIC ARMA_DONT_USE_WRAPPER)
#target_compile_definitions(smoldock_object PUBLIC BOOST_LOG_DYN_LINK)



add_custom_command(
        OUTPUT partially_linked_smoldock.o
        COMMAND ld -Ur $<TARGET_OBJECTS:smoldock_object>   -o partially_linked_smoldock.o -Bstatic -L${VC_ROOT}/lib/ -lVc -L${RDKIT_ROOT}/lib -L ${PATH_TO_COMPILER_INTERNAL_LIBS_1} -L /usr/lib -L ${PATH_TO_COMPILER_INTERNAL_LIBS_2} -L ${PATH_TO_COMPILER_INTERNAL_LIBS_3}  -lgfortran -llapack -lblas -lgfortran -lquadmath  -lgfortran   -lRDKitStatic -lunwind -lRDKitForceFieldHelpers_static -lRDKitForceField_static -lRDKitStatic -L${Boost_LIBRARY_DIRS} -lboost_log -lboost_log_setup  -lboost_system -lboost_serialization -lboost_log -lboost_system -lboost_filesystem -lboost_date_time -lboost_thread -lboost_regex -lboost_mpi
        DEPENDS $<TARGET_OBJECTS:smoldock_object>
        COMMAND_EXPAND_LISTS
)



ADD_LIBRARY(libsmoldock STATIC partially_linked_smoldock.o)
set_target_properties(libsmoldock PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(libsmoldock PROPERTIES COMPILE_FLAGS "")
set_target_properties(libsmoldock PROPERTIES LINK_FLAGS " -static -v  ")
set_target_properties(libsmoldock PROPERTIES PREFIX "lib")
set_target_properties(libsmoldock PROPERTIES OUTPUT_NAME "smoldock")
target_compile_definitions(libsmoldock PUBLIC ARMA_DONT_USE_WRAPPER)
add_dependencies(libsmoldock smoldock_object)

ADD_LIBRARY(libsmoldockPIC SHARED ${SOURCE_LIST} Engines/ScoringFunctions/VinaExtendedCommon.h Utilities/LogUtils.h)
set_target_properties(libsmoldockPIC PROPERTIES COMPILE_FLAGS " -fPIC ")
set_target_properties(libsmoldockPIC PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(libsmoldockPIC PROPERTIES LINK_FLAGS "  -shared  -fPIC   ")
set_target_properties(libsmoldockPIC PROPERTIES PREFIX "lib")
set_target_properties(libsmoldockPIC PROPERTIES OUTPUT_NAME "smoldockPIC")
target_compile_definitions(libsmoldockPIC PUBLIC BOOST_LOG_DYN_LINK)
target_compile_definitions(libsmoldockPIC PUBLIC ARMA_DONT_USE_WRAPPER)
target_link_libraries(libsmoldockPIC -Wl,-Bdynamic -Wl,-lboost_regex  -Wl,-lboost_filesystem)
target_link_libraries(libsmoldockPIC -Wl,-Bdynamic  -Wl,-lboost_log   -Wl,-lboost_system  -Wl,-lboost_thread  )
target_link_libraries(libsmoldockPIC -Wl,-Bdynamic -Wl,-lboost_numpy -Wl,-lboost_python-3.6-mt  -Wl,-lpython3.6m)
TARGET_LINK_LIBRARIES(libsmoldockPIC -Wl,-Bstatic -Wl,-lRDKitStatic -Wl,-lRDKitForceFieldHelpers_static -Wl,-lRDKitForceField_static)
target_link_libraries(libsmoldockPIC -Wl,-Bstatic -Wl,-lboost_date_time  -Wl,-lboost_date_time)
target_link_libraries(libsmoldockPIC -Wl,-Bstatic -Wl,-lVc)
target_link_libraries(libsmoldockPIC -Wl,-Bstatic -Wl,-llapack -Wl,-lblas -Wl,-l:libquadmath.a)
target_link_libraries(libsmoldockPIC -Wl,-Bdynamic -Wl,-lgfortran)
target_link_libraries(libsmoldockPIC -Wl,-Bdynamic -Wl,-L/usr/lib/gcc/x86_64-pc-linux-gnu/7.3.0/ -Wl,-lgcc)
target_link_libraries(libsmoldockPIC -Wl,-Bdynamic -Wl,-lpthread -Wl,-ldl)

#add_compile_definitions(BOOST_LOG_DYN_LINK BOOST_TEST_DYN_LINK)


#target_compile_definitions(smoldock PUBLIC BOOST_LOG_DYN_LINK)


add_executable(smoldock Frontends/main.cpp)
set_target_properties(smoldock PROPERTIES LINK_FLAGS " -Wl,--no-undefined -Wl,--as-needed -static  -Wl,-allow-multiple-definition   ")
TARGET_LINK_LIBRARIES(smoldock libsmoldock)
TARGET_LINK_LIBRARIES(smoldock pthread)
target_link_libraries(smoldock boost_log bprinter tbb_static)

add_executable(smoldock_dyn Frontends/main.cpp)
set_target_properties(smoldock_dyn PROPERTIES LINK_FLAGS " -Wl,--no-undefined -Wl,--as-needed    ")
TARGET_LINK_LIBRARIES(smoldock_dyn libsmoldock)
TARGET_LINK_LIBRARIES(smoldock_dyn pthread)
target_link_libraries(smoldock_dyn -Wl,-Bstatic boost_log bprinter tbb_static  -Wl,-Bdynamic)

if(SMOLDOCK_INTEL_VTUNE_EXPORT_SYMBOL)
    target_link_options(smoldock PRIVATE "-u_init" "-usetenv" "-ugetenv" "-u__errno_location" "-upthread_create" "-upthread_getattr_np")
    target_link_options(smoldock PRIVATE "-upthread_self" "-upthread_getspecific" "-upthread_setspecific")
    target_link_options(smoldock PRIVATE "-upthread_getattr_np" "-upthread_attr_destroy" "-upthread_attr_setstack" "-upthread_attr_getstack")
    target_link_options(smoldock PRIVATE "-upthread_attr_getstacksize" "-upthread_attr_setstacksize" "-upthread_cancel" "-u_pthread_cleanup_push")
    target_link_options(smoldock PRIVATE "-u_pthread_cleanup_pop" "-upthread_attr_setstacksize" "-upthread_cancel" "-u_pthread_cleanup_push")
    target_link_options(smoldock PRIVATE "-upthread_mutex_lock" "-upthread_mutex_trylock" "-upthread_spin_lock" "-upthread_spin_trylock")
endif()


add_executable(goldenRecordGenerator Frontends/goldenRecordGenerator.cpp)
set_target_properties(goldenRecordGenerator PROPERTIES LINK_FLAGS " -Wl,--no-undefined -Wl,--as-needed -static  -Wl,-allow-multiple-definition   ")
TARGET_LINK_LIBRARIES(goldenRecordGenerator libsmoldock)
TARGET_LINK_LIBRARIES(goldenRecordGenerator pthread)
target_link_libraries(goldenRecordGenerator boost_log bprinter tbb_static)


add_executable(mpicalibrator Frontends/mpicalibrator.cpp ${MPI_SOURCE_LIST})
set_target_properties(mpicalibrator PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS} ")
set_target_properties(mpicalibrator PROPERTIES LINK_FLAGS " -Wl,--as-needed -Wl,-allow-multiple-definition ")
TARGET_LINK_LIBRARIES(mpicalibrator -L${PATH_TO_COMPILER_INTERNAL_LIBS_1} -L${PATH_TO_COMPILER_INTERNAL_LIBS_2} -L ${PATH_TO_COMPILER_INTERNAL_LIBS_3})
TARGET_LINK_LIBRARIES(mpicalibrator -Wl,-Bstatic libsmoldock )
target_link_libraries(mpicalibrator -Wl,-Bstatic boost_log boost_mpi boost_serialization bprinter tbb_static  )
TARGET_LINK_LIBRARIES(mpicalibrator -Wl,-Bdynamic pthread  ${MPI_LINK_FLAGS})
target_link_libraries(mpicalibrator -Wl,-Bstatic -Wl,-llapack -Wl,-lblas -Wl,-l:libquadmath.a)
target_link_libraries(mpicalibrator -Wl,-Bdynamic dl -Wl,-lgfortran)
target_link_libraries(mpicalibrator -Wl,-Bstatic Vc -Wl,-Bdynamic)




add_executable(vina_refiner Frontends/vina_refiner.cpp)
set_target_properties(vina_refiner PROPERTIES LINK_FLAGS " -Wl,--no-undefined -Wl,--as-needed  -static   ")
TARGET_LINK_LIBRARIES(vina_refiner libsmoldock)
TARGET_LINK_LIBRARIES(vina_refiner pthread)
target_link_libraries(vina_refiner boost_log boost_filesystem boost_program_options)

add_executable(smoldockCLI Frontends/cli_interface.cpp)
set_target_properties(smoldockCLI PROPERTIES LINK_FLAGS " -Wl,--no-undefined -Wl,--as-needed  -static   ")
TARGET_LINK_LIBRARIES(smoldockCLI libsmoldock)
TARGET_LINK_LIBRARIES(smoldockCLI pthread)
target_link_libraries(smoldockCLI boost_log boost_filesystem boost_program_options)




add_library(PySmolDock MODULE Frontends/SmoldockPythonBindings.cpp ${PY_EXPORT_SOURCE_LIST})
set_target_properties(PySmolDock PROPERTIES COMPILE_FLAGS " -fPIC ")
set_target_properties(PySmolDock PROPERTIES LINK_FLAGS " -shared -Wl,-nostdlib -Wl,-L./ -Wl,-lsmoldockPIC -Wl,-Bshareable")
set_target_properties(PySmolDock PROPERTIES PREFIX "")
set_target_properties(PySmolDock PROPERTIES OUTPUT_NAME "PySmolDock")
target_compile_definitions(PySmolDock PUBLIC BOOST_LOG_DYN_LINK)
target_compile_definitions(PySmolDock PUBLIC ARMA_DONT_USE_WRAPPER)
target_link_libraries(PySmolDock )
add_dependencies(PySmolDock libsmoldockPIC)
### Remove this for prod build

IF(CMAKE_BUILD_TYPE MATCHES Debug)
    target_compile_definitions(libsmoldock PUBLIC SMOLDOCK_VERBOSE_DEBUG)
    target_compile_definitions(smoldock PUBLIC SMOLDOCK_VERBOSE_DEBUG)
    target_compile_definitions(vina_refiner PUBLIC SMOLDOCK_VERBOSE_DEBUG)
    target_compile_definitions(PySmolDock PUBLIC SMOLDOCK_VERBOSE_DEBUG)
    target_compile_definitions(smoldock_object PUBLIC SMOLDOCK_VERBOSE_DEBUG)
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
TARGET_LINK_LIBRARIES(complete_test_suite libsmoldock)
TARGET_LINK_LIBRARIES(complete_test_suite pthread)
target_link_libraries(complete_test_suite boost_log)
target_link_libraries(complete_test_suite bprinter tbb_static)
target_link_libraries(complete_test_suite ${Boost_LIBRARIES} -Wl,-Bdynamic boost_unit_test_framework )
#target_link_libraries(complete_test_suite Vc)


