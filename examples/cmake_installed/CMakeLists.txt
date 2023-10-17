cmake_minimum_required(VERSION 3.18.0)
project(pressio-cmake-example)

find_package(pressio REQUIRED)

add_executable(pressio-cmake-example main.cpp)
target_link_libraries(pressio-cmake-example PUBLIC pressio)
