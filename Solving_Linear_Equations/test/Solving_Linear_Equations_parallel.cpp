#include "Solving_Linear_Equations.h"
#include <iostream>

#define FIRST_THREAD 0

void Solving_Linear_Equations_parallel::proximity_function()
{
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int count_for_process = ceil((double)N / size);
	int ibeg = count_for_process * rank;
	int iend = count_for_process * (rank + 1) > N ? N : count_for_process * (rank + 1);

	std::vector<double> x_process(iend - ibeg);
	for (int i = ibeg, j = 0; i < iend; i++, j++) // ллллл
	{
		x_process[j] = multiply_row_by_column(A[i], x);
	}
	for (int i = ibeg, j = 0; i < iend; i++, j++)
	{
		x_process[j] = x[i] - ti * (x_process[j] - b[i]);
	}

	if (rank == FIRST_THREAD)
	{
		for (int i = 0; i < count_for_process; i++)
		{
			x[i] = x_process[i];
		}
		for(int i = 1; i < size; i++)
		{
			MPI_Recv(&x[i * count_for_process], count_for_process, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		for (int i = 1; i < size; i++)
		{
			MPI_Ssend(&x[0], x.size(), MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
	}
	else
	{
		MPI_Ssend(&x_process[0], x_process.size(), MPI_DOUBLE, FIRST_THREAD, 0, MPI_COMM_WORLD);
		MPI_Recv(&x[0], x.size(), MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
}

bool Solving_Linear_Equations_parallel::accuracy_check(double epsilon) const
{
	std::vector<double> result(N);
	for (int i = 0; i < N; i++)
	{
		result[i] = multiply_row_by_column(A[i], x);
	}

	for (int i = 0; i < N; i++)
	{
		result[i] = (result[i] - b[i]);
	}

	double norm_numerator = find_norm(result);
	return norm_numerator / norm_denominator < epsilon ? true : false;
}
