add_executable(test test.cpp)
target_link_libraries(test gtest gtest_main fl_rmq_lib)

add_executable(test_bp test_bp.cpp)
target_link_libraries(test_bp gtest gtest_main fl_rmq_lib)

# Add the test executable to CTest
#add_test(NAME test COMMAND test)
