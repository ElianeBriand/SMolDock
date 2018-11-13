project(unit_tests)
enable_testing()


find_package(Boost REQUIRED COMPONENTS unit_test_framework)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})

add_executable(main_test_suite main_test_suite.cpp $<TARGET_OBJECTS:smoldock_object>)

target_link_libraries(main_test_suite ${Boost_LIBRARIES})