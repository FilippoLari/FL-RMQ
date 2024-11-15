cmake_minimum_required(VERSION 3.12)
project(fl_rmq)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -g -O3")

set(CMAKE_CXX_STANDARD 20)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if (COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif ()

include_directories(include)

add_library(fl_rmq_lib INTERFACE)
target_include_directories(fl_rmq_lib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# add libsais
set(LIBSAIS_DIR ${CMAKE_SOURCE_DIR}/libsais)
set(LIBSAIS_LIB ${LIBSAIS_DIR}/liblibsais.a)

add_library(libsais STATIC IMPORTED)
set_target_properties(libsais PROPERTIES IMPORTED_LOCATION ${LIBSAIS_LIB})

target_include_directories(libsais INTERFACE ${LIBSAIS_DIR}/include)

add_subdirectory(tests)

add_subdirectory(experiments)

add_subdirectory(googletest)
#enable_testing()