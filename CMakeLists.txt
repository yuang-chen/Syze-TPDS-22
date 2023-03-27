
cmake_minimum_required(VERSION 3.22)

project(syze VERSION 0.1.0)

find_package(fmt CONFIG REQUIRED)
find_package(OpenMP REQUIRED)

add_library(common INTERFACE)
set_target_properties(common PROPERTIES CXX_EXTENSIONS OFF)

target_compile_features(common INTERFACE cxx_std_20)
target_compile_options(common INTERFACE -fopenmp)

add_compile_options("$<$<CONFIG:Debug>:-DDEBUG>")

# include_directories()
target_link_libraries(common INTERFACE fmt::fmt)
target_link_libraries(common INTERFACE OpenMP::OpenMP_CXX)

add_executable(pr app/pr.cpp)
add_executable(cf app/cf.cpp)
add_executable(ec app/ec.cpp)
add_executable(in app/indegree.cpp)

# target_compile_options(test PRIVATE -std=c++20 -ftree-vectorize -mavx512f -mavx512cd -ffast-math -march=skylake-avx512)
target_link_libraries(pr PRIVATE common)
target_link_libraries(cf PRIVATE common)
target_link_libraries(ec PRIVATE common)
target_link_libraries(in PRIVATE common)
