cmake_minimum_required(VERSION 3.14)

set(CMAKE_CXX_STANDARD 17)

project(Solving_Linear_Equations)

set(exe_name ${PROJECT_NAME})

#add_library(corelib OBJECT word_counter.cpp)

add_executable(${exe_name} main.c)
#target_link_libraries(${exe_name} corelib)
