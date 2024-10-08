# test/CMakeLists.txt

file(GLOB TEST_SOURCES *.cpp)
add_executable(run_tests ${TEST_SOURCES})

# Link with Google Test
target_link_libraries(run_tests gtest gtest_main)


add_test(NAME run_tests COMMAND run_tests)
