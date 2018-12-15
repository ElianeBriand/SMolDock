cmake_minimum_required(VERSION 3.12)
project(smoldock)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})

SET (CMAKE_CXX_FLAGS                "-Wall")
SET (CMAKE_CXX_COMPILER             "clang++")
SET (CMAKE_CXX_FLAGS_DEBUG          "-g")
SET (CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
SET (CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

SET (CMAKE_AR      "llvm-ar")
SET (CMAKE_LINKER  "llvm-ld")
SET (CMAKE_NM      "llvm-nm")
SET (CMAKE_OBJDUMP "llvm-objdump")
SET (CMAKE_RANLIB  "llvm-ranlib")



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
        Engines/ScoringFunctions/BasicScoringFunction.cpp
        Engines/ScoringFunctions/BasicScoringFunction.h
        Engines/Internals/iProtein.h Engines/Internals/iTransform.h)

ADD_LIBRARY(smoldock_object OBJECT ${SOURCE_LIST})

### Remove this for prod build
add_compile_definitions(SMOLDOCK_VERBOSE_DEBUG)

add_executable(smoldock Frontends/main.cpp $<TARGET_OBJECTS:smoldock_object>)
set_target_properties(smoldock PROPERTIES LINK_FLAGS " -static  ")


TARGET_LINK_LIBRARIES(smoldock ${RDKIT_LINKAGE_LIST})


target_link_libraries(smoldock gsl gslcblas boost_log boost_log_setup boost_system boost_log boost_system boost_filesystem boost_date_time boost_thread boost_regex pthread)
target_link_libraries(smoldock Vc)


###############################################################################################################
############                                  UNIT TESTING
###############################################################################################################
enable_testing()

set(TEST_SOURCE_LIST
        UnitTests/main_test_suite.cpp
        UnitTests/UnitTestHelper.cpp
        UnitTests/UnitTestHelper.h)


add_executable(main_test_suite ${TEST_SOURCE_LIST} $<TARGET_OBJECTS:smoldock_object>)

set_target_properties(main_test_suite PROPERTIES LINK_FLAGS " -static  ")



TARGET_LINK_LIBRARIES(main_test_suite ${RDKIT_LINKAGE_LIST})

target_link_libraries(main_test_suite gsl gslcblas boost_log boost_log_setup  boost_system boost_system boost_filesystem boost_date_time boost_thread boost_regex )
target_link_libraries(main_test_suite Vc)

target_link_libraries(main_test_suite ${Boost_LIBRARIES} boost_unit_test_framework)
