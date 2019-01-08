cmake_minimum_required(VERSION 3.12)
project(smoldock)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost REQUIRED)

set(PY_VERSION 3.6)

find_package(PythonLibs ${PY_VERSION} REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})

SET (CMAKE_CXX_FLAGS                "-Wall")
SET (CMAKE_CXX_COMPILER             "clang++")
SET (CMAKE_CXX_FLAGS_DEBUG          "-g")
SET (CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
SET (CMAKE_CXX_FLAGS_RELEASE        "-O3 -march=native -DNDEBUG")
SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -march=native -g ")



######## RDKIT SETUP ######
# Adjust this according to where your RDKit build is :
set(RDKIT_ROOT /home/eliane/Projects/RDKit_Git/rdkit-Release_2018_09_1/)

include_directories(${RDKIT_ROOT}/Code/)

LINK_DIRECTORIES(${RDKIT_ROOT}/lib/staticlib)
LINK_DIRECTORIES(${RDKIT_ROOT}/lib)

set(RDKIT_LINKAGE_LIST RDKitStatic pthread RDKitStatic RDKitForceField_static RDKitForceFieldHelpers_static)
######## END OF RDKIT ######


######## VC SETUP ######
# Adjust this according to where your Vc build is :
set(VC_ROOT /home/eliane/Projects/VcInstall/)

include_directories(${VC_ROOT}/include/)

LINK_DIRECTORIES(${VC_ROOT}/lib)
######## END OF VC  ######

include_directories(Dependencies)
include_directories(.)


include_directories(${PYTHON_INCLUDE_DIRS})

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
        Engines/ConformerRigidDockingEngine.cpp
        Engines/ConformerRigidDockingEngine.h
        Engines/Internals/iConformer.h
        Utilities/TimingsLog.h
        Engines/ScoringFunctions/VinaLikeScoringFunction.cpp
        Engines/ScoringFunctions/VinaLikeScoringFunction.h
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
        Structures/Common/VariantFlagAssignation.cpp
        Structures/Common/VariantFlagAssignation.h
        Structures/Common/PDBResidueVariantTable.cpp
        Structures/Common/PDBResidueVariantTable.h
        Structures/InputPostProcessors/VinaCompatibilityPostProcessor.cpp
        Structures/InputPostProcessors/VinaCompatibilityPostProcessor.h
        Structures/InputPostProcessors/InputPostProcessorInterface.h
        Engines/ScoringFunctions/ScoringFunctionInterface.h
        Engines/LocalOptimizers/L_BFGS.cpp
        Engines/LocalOptimizers/L_BFGS.h
        Engines/LocalOptimizers/OptimizerInterface.h
        Engines/Internals/InternalsUtilityFunctions.h
        Engines/LocalOptimizers/LocalOptimizerUtilityFunctions.h
        Engines/GlobalHeuristics/RandomRestart.cpp
        Engines/GlobalHeuristics/RandomRestart.h
        Engines/GlobalHeuristics/HeuristicInterface.h Engines/ScoringFunctions/ScoringFunctionFactory.cpp Engines/ScoringFunctions/ScoringFunctionFactory.h)

set(PY_EXPORT_SOURCE_LIST Frontends/Python/PyStructures.cpp Frontends/Python/PyStructures.h Frontends/Python/PySTLWrapper.cpp Frontends/Python/PySTLWrapper.h Frontends/Python/PyUtilities.cpp Frontends/Python/PyUtilities.h Frontends/Python/PyScoringFunctions.cpp Frontends/Python/PyScoringFunctions.h)

ADD_LIBRARY(smoldock_object OBJECT ${SOURCE_LIST})
set_target_properties(smoldock_object PROPERTIES LINK_FLAGS " -static  ")
target_compile_definitions(smoldock_object PUBLIC ARMA_DONT_USE_WRAPPER)

ADD_LIBRARY(libsmoldock STATIC ${SOURCE_LIST} )
set_target_properties(libsmoldock PROPERTIES COMPILE_FLAGS " -fPIC ")
set_target_properties(libsmoldock PROPERTIES LINK_FLAGS " -static -v -Wl,--no-undefined -Wl,--as-needed ")
set_target_properties(libsmoldock PROPERTIES PREFIX "lib")
set_target_properties(libsmoldock PROPERTIES OUTPUT_NAME "smoldock")
TARGET_LINK_LIBRARIES(libsmoldock ${RDKIT_LINKAGE_LIST})
target_link_libraries(libsmoldock gsl gslcblas boost_log boost_log_setup boost_system boost_log boost_system boost_filesystem boost_date_time boost_thread boost_regex pthread)
target_link_libraries(libsmoldock Vc)
target_link_libraries(libsmoldock blas lapack gfortran quadmath)
target_compile_definitions(libsmoldock PUBLIC BOOST_LOG_DYN_LINK)
target_compile_definitions(smoldock_object PUBLIC ARMA_DONT_USE_WRAPPER)

### Remove this for prod build
add_compile_definitions(SMOLDOCK_VERBOSE_DEBUG)

# You need this if you do NOT generate a static executable
#add_compile_definitions(BOOST_LOG_DYN_LINK BOOST_TEST_DYN_LINK)

add_executable(smoldock Frontends/main.cpp $<TARGET_OBJECTS:smoldock_object>)

# Comment this out if you do NOT want to build a static executable
set_target_properties(smoldock PROPERTIES LINK_FLAGS " -Wl,--no-undefined  -Wl,--as-needed  -static  ")


TARGET_LINK_LIBRARIES(smoldock ${RDKIT_LINKAGE_LIST})


target_link_libraries(smoldock gsl gslcblas boost_log boost_log_setup boost_system boost_log boost_system boost_filesystem boost_date_time boost_thread boost_regex pthread)
target_link_libraries(smoldock Vc)
target_link_libraries(smoldock blas lapack gfortran quadmath)



# For Intel VTune/Profiler tools, you need this if you generate a static executable
#target_link_options(smoldock PRIVATE "-u_init" "-usetenv" "-ugetenv" "-u__errno_location" "-upthread_create" "-upthread_getattr_np")
#target_link_options(smoldock PRIVATE "-upthread_self" "-upthread_getspecific" "-upthread_setspecific")
#target_link_options(smoldock PRIVATE "-upthread_getattr_np" "-upthread_attr_destroy" "-upthread_attr_setstack" "-upthread_attr_getstack")
#target_link_options(smoldock PRIVATE "-upthread_attr_getstacksize" "-upthread_attr_setstacksize" "-upthread_cancel" "-u_pthread_cleanup_push")
#target_link_options(smoldock PRIVATE "-u_pthread_cleanup_pop" "-upthread_attr_setstacksize" "-upthread_cancel" "-u_pthread_cleanup_push")
#target_link_options(smoldock PRIVATE "-upthread_mutex_lock" "-upthread_mutex_trylock" "-upthread_spin_lock" "-upthread_spin_trylock")

add_library(PySmolDock MODULE  Frontends/SmoldockPythonBindings.cpp ${PY_EXPORT_SOURCE_LIST} )
set_target_properties(PySmolDock PROPERTIES COMPILE_FLAGS " -fPIC ")
set_target_properties(PySmolDock PROPERTIES LINK_FLAGS " -Wl,--no-undefined  -Wl,--as-needed ")
set_target_properties(PySmolDock PROPERTIES PREFIX "")
set_target_properties(PySmolDock PROPERTIES OUTPUT_NAME "PySmolDock")

# link
target_link_libraries(PySmolDock libsmoldock boost_numpy boost_python-3.6-mt python3.6m boost_log_setup boost_log boost_date_time)



###############################################################################################################
############                                  UNIT TESTING
###############################################################################################################
enable_testing()

set(TEST_SOURCE_LIST
        UnitTests/main_test_suite.cpp
        UnitTests/engine_test_suite.cpp)


add_executable(complete_test_suite ${TEST_SOURCE_LIST} $<TARGET_OBJECTS:smoldock_object>)

# Comment this out if you do NOT want to build a static executable
set_target_properties(complete_test_suite PROPERTIES LINK_FLAGS " -static  ")



TARGET_LINK_LIBRARIES(complete_test_suite ${RDKIT_LINKAGE_LIST})

target_link_libraries(complete_test_suite gsl gslcblas boost_log boost_log_setup  boost_system boost_system boost_filesystem boost_date_time boost_thread boost_regex )
target_link_libraries(complete_test_suite Vc)
target_link_libraries(complete_test_suite blas lapack gfortran quadmath)
target_link_libraries(complete_test_suite ${Boost_LIBRARIES} boost_unit_test_framework)
