cmake_minimum_required(VERSION 3.5)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(pybind11_DIR "../venv/Lib/site-packages/pybind11/share/cmake/pybind11")
project(model VERSION 1.0)
find_package(pybind11 REQUIRED)
pybind11_add_module(model SHARED binds.cpp channel.cpp shadowing.cpp)
target_include_directories(model PUBLIC "${CMAKE_SOURCE_DIR}")

