#include<mpi.h>
#include<vector>
#include<iostream>
#include <chrono>


#include "Solving_Linear_Equations.h"

int main(int argc, char* argv[])
{
	auto start_time = std::chrono::high_resolution_clock::now();
	MPI_Init(&argc, &argv);
	int N = 5000;
	Matrix A(N);
	std::vector<double> x(N, 0.0);
	std::vector<double> b(N, N + 1);
	std::vector<double> res(N);
	Solving_Linear_Equations_parallel test(A, x, b);
	res = test.execute(0.00001);
	auto end_time = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
	std::cout << "Time passed: " << duration.count() << " micsec" << std::endl;
	MPI_Finalize();
}