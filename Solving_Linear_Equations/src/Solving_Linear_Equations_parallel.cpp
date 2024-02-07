#include "Solving_Linear_Equations.h"

#define FIRST_THREAD 0

void Solving_Linear_Equations_parallel::proximity_function()
{
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int count_for_process = ceil((double)N / size);
	int ibeg = count_for_process * rank;
	int iend = count_for_process * (rank + 1) > N ? N : count_for_process * (rank + 1);

	std::vector<double> x_process(count_for_process);
	for (int i = ibeg, j = 0; i < iend; i++, j++)
	{
		x_process[j] = multiply_row_by_column(A[i], x);
	}
	for (int i = ibeg, j = 0; i < iend; i++, j++)
	{
		x_process[j] = x[i] - ti_plus * (x_process[j] - b[i]);
	}

	if (rank == FIRST_THREAD)
	{
		std::vector<double> new_x(N);
		new_x = x_process;
		if (size > 1)
		{
			MPI_Recv(&x_process + ibeg, x_process.size(), MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		x = new_x;
	}
	else
	{
		MPI_Send(&x_process, x_process.size(), MPI_DOUBLE, FIRST_THREAD, 0, MPI_COMM_WORLD);
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
