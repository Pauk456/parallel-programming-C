#!/bin/bash
mpic++ main.cpp Solving_Linear_Equations_virtual.cpp Solving_Linear_Equations_usual.cpp Solving_Linear_Equations_parallel_first.cpp Solving_Linear_Equations_parallel_second.cpp Matrix.cpp -o main
num_processes=${1:-1}
mpiexec -n $num_processes ./main
rm main
