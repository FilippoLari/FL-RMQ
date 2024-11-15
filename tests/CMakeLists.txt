add_executable(test test.cpp)

target_link_libraries(test gtest gtest_main sdsl fl_rmq_lib)

# Add the test executable to CTest
#add_test(NAME test COMMAND test)
