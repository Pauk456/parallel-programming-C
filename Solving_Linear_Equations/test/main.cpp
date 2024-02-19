#include<mpi.h>
#include<vector>
#include<iostream>
#include<chrono>
#include "Solving_Linear_Equations.h"

int main(int argc, char* argv[])
{
	int N = 1000;

	Matrix A(N);
	std::vector<double> x(N, 2.0);
	std::vector<double> b(N, N + 1);
	try {
		auto start_time = std::chrono::high_resolution_clock::now();

		Solving_Linear_Equations_parallel_second solver(A, x, b, argc, argv);
		solver.execute(0.0001);

		auto end_time = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
		std::cout << "Time passed: " << duration.count() << " micsec" << std::endl;

		solver.print_result();
	}
	catch (const std::exception& e)
	{
		std::cout << e.what() << std::endl;
	}
}
