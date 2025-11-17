cmake_minimum_required(VERSION 3.16)
project(division_period)

set(CMAKE_CXX_STANDARD 20)

add_executable(division_period normal.cpp)

target_link_libraries(division_period pthread)

add_executable(simple simple.cpp)