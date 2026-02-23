cmake_minimum_required(VERSION 3.16)
project(division_period LANGUAGES CXX)

# ---- C++ 标准 ----
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---- 编译警告 ----
add_compile_options(-Wall -Wextra -Wpedantic)

# ---- 线程库（跨平台方式） ----
find_package(Threads REQUIRED)

# ---- 可执行文件：数论版 ----
add_executable(division_period normal.cpp)
target_link_libraries(division_period PRIVATE Threads::Threads)

# ---- 可执行文件：简洁版 ----
add_executable(simple simple.cpp)
