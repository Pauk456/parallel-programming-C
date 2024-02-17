cmake_minimum_required(VERSION 3.14)

set(CMAKE_CXX_STANDARD 17)

project(Solving_Linear_Equations)

set(exe_name ${PROJECT_NAME})

set(MPI_ROOT "C://Program Files (x86)//Microsoft SDKs//MPI")

include_directories(${MPI_ROOT}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

link_directories(${MPI_ROOT}/Lib/x64)

add_executable(
	${exe_name} main.cpp 
	src/Solving_Linear_Equations_usual.cpp 
	src/Matrix.cpp 
	src/Solving_Linear_Equations_virtual.cpp
	src/Solving_Linear_Equations_parallel_first.cpp
	src/Solving_Linear_Equations_parallel_second.cpp
)
target_link_libraries(${exe_name} msmpi.lib)
