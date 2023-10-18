cmake_minimum_required(VERSION 3.18.0)
project(Example)

add_subdirectory(${Example_SOURCE_DIR}/../.. ${Example_BINARY_DIR}/pressio)

add_executable(pressio-cmake-example main.cpp)
target_link_libraries(pressio-cmake-example PUBLIC pressio)
target_compile_definitions(pressio-cmake-example PUBLIC PRESSIO_ENABLE_TPL_EIGEN)
