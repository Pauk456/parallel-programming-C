#include<mpi.h>
#include<vector>
#include<iostream>
#include<chrono>


#include "Solving_Linear_Equations.h"

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	auto start_time = std::chrono::high_resolution_clock::now();
	int N = 5000;
	Matrix A(N);
	std::vector<double> x(N, 2.0);
	std::vector<double> b(N, N + 1);
	std::vector<double> res(N);
	int rank;
	Solving_Linear_Equations_parallel test(A, x, b);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	res = test.execute(0.0001);
	auto end_time = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
	if (rank == 0)
	{
		std::cout << "Time passed: " << duration.count() << " micsec" << "rank:" << rank << std::endl;
		//for (int i = 0; i < N; i++)
		//{
		//	std::cout << res[i] << ' ';
		//}
		//std::cout << std::endl;
	}
}
