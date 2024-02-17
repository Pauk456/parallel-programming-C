#include "Solving_Linear_Equations.h"

#define FIRST_THREAD 0

Solving_Linear_Equations_parallel_second::Solving_Linear_Equations_parallel_second(Matrix A, std::vector<double> x, std::vector<double> b)
: Solving_Linear_Equations_virtual(A)
{
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	count_for_process = ceil((double)N / size);
	ibeg = count_for_process * rank;
	iend = count_for_process * (rank + 1) > N ? N : count_for_process * (rank + 1);
	count_for_process = iend - ibeg;

	x.resize(count_for_process);
	b.resize(count_for_process);
	for (int i = ibeg; i < iend; i++)
	{
		this->x[i] = x[i];
		this->b[i] = b[i];
	}


};

double Solving_Linear_Equations_parallel_second::multiply_row_by_column(const std::vector<double>& row, const std::vector<double>& column, int offset, int count) const
{
	double result = 0.0;
	for (int i = offset; i < count; ++i) {
		result += row[i] * column[i];
	}
	return result;
}

void Solving_Linear_Equations_parallel_second::proximity_function()
{
	std::vector<double> x_process(count_for_process, 0.0);

	for (int k = 0; k < count_for_process; k++)
	{
		int block = (rank + k) % size;
		
		for (int i = ibeg, j = 0; i < iend; i++, j++) // �����
		{
			x_process[j] += multiply_row_by_column(A[i], x, block * count_for_process,count_for_process);
		}

	}
	


	for (int i = ibeg, j = 0; i < iend; i++, j++)
	{
		x_process[j] = x[i] - ti * (x_process[j] - b[i]);
	}


	int process_collector = 0;
	if (rank == FIRST_THREAD)
	{
		for (int i = 0; i < count_for_process; i++)
		{
			x[i] = x_process[i];
		}
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
		MPI_Ssend(&x_process[0], x_process.size(), MPI_DOUBLE, process_collector, 0, MPI_COMM_WORLD);
		MPI_Recv(&x[0], x.size(), MPI_DOUBLE, process_collector, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
}

bool Solving_Linear_Equations_parallel_second::accuracy_check(double epsilon) const
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
	return norm_numerator / norm_denominator < epsilon * epsilon ? true : false;
}