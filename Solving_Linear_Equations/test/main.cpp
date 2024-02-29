#include<mpi.h>
#include<vector>
#include<iostream>
#include<chrono>
#include "Solving_Linear_Equations.h"

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int N = 1000;

	Matrix A(N);
	std::vector<double> x(N, 0);
	std::vector<double> b(N, N + 1);
	try {
		auto start_time = std::chrono::high_resolution_clock::now();

		Solving_Linear_Equations_parallel_second solver(A, x, b);
		solver.execute(0.0001);

		auto end_time = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
		if (rank == 0)
		{
			std::cout << "Time passed: " << duration.count() << " micsec" << std::endl;
		}
		solver.print_result();
	}
	catch (const std::exception& e)
	{
		std::cout << e.what() << std::endl;
	}
	
	MPI_Finalize();
}
