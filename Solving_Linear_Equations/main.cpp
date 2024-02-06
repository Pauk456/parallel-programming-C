#include<mpi.h>
#include<vector>
#include<iostream>

#include "Solving_Linear_Equations.h"

int main(int argc, char* argv[])
{
	int N = 50;
	Matrix A(N);
	std::vector<double> x(N, 0.0);
	std::vector<double> b(N, N + 1);
	std::vector<double> res(N);
	Solving_Linear_Equations test(A, x, b);
	res = test.execute(0.00001);
	for (int i = 0; i < N; i++)
	{
		std::cout << res[i] << ' ';
	}
}