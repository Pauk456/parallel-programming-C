#include "Solving_Linear_Equations.h"
#include<mpi.h>

Solving_Linear_Equations_parallel_first::Solving_Linear_Equations_parallel_first(Matrix A, std::vector<double> x, std::vector<double> b)
	: Solving_Linear_Equations_virtual(A, x, b)
{
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (N % size != 0)
	{
		throw std::runtime_error("Error: N % size is not 0");
	}

	count_for_process = ceil((double)N / size);
	ibeg = count_for_process * rank;
	iend = count_for_process * (rank + 1) > N ? N : count_for_process * (rank + 1);
	count_for_process = iend - ibeg;

	result_part.resize(count_for_process);
	result.resize(N);
};

void Solving_Linear_Equations_parallel_first::proximity_function()
{ 
	for (int i = ibeg; i < iend; i++)
	{
		x[i] = x[i] - ti * (multiply_row_by_column(A[i], x) - b[i]);
	}

	//MPI_Allgather(&x[ibeg], count_for_process, MPI_DOUBLE, &x[0], count_for_process, MPI_DOUBLE, MPI_ANY_TAG);
	int process_collector = 0;
	if (rank == process_collector)
	{
		for (int i = 1; i < size; i++)
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
		MPI_Ssend(&x[ibeg], count_for_process, MPI_DOUBLE, process_collector, 0, MPI_COMM_WORLD);
		MPI_Recv(&x[0], x.size(), MPI_DOUBLE, process_collector, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
}

bool Solving_Linear_Equations_parallel_first::accuracy_check(double epsilon)
{
	for (int i = ibeg, j = 0; i < iend; i++, j++)
	{
		result_part[j] = multiply_row_by_column(A[i], x) - b[i];
	}

	MPI_Allgather(result_part.data(), count_for_process, MPI_DOUBLE, result.data(), count_for_process, MPI_DOUBLE, MPI_COMM_WORLD);

	double norm_numerator = find_norm(result);
	return norm_numerator / norm_denominator < epsilon;
}

void Solving_Linear_Equations_parallel_first::print_result()
{
	if (rank == 0)
	{
		std::cout << "Count of elements X = " << x.size() << std::endl;
		for (int i = 0; i < x.size(); i++)
		{
			std::cout << x[i] << ' ';
		}
		std::cout << std::endl;
	}
}
