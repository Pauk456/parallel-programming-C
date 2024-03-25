cmake_minimum_required(VERSION 3.14)

set(CMAKE_CXX_STANDARD 17)

project(Linear_Equations_OpenMP)

set(exe_name ${PROJECT_NAME})

add_executable(
	${exe_name} First.c Second.c Usual.c  "a.c" "Second2.c" "Second1.cpp")
