add_executable(build_uniform build_uniform.cpp)
target_link_libraries(build_uniform fl_rmq_lib)

add_executable(bar_plot_lcp bar_plot_lcp.cpp)
target_link_libraries(bar_plot_lcp libsais fl_rmq_lib)
