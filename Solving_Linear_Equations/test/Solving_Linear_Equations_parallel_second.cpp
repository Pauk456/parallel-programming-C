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

	destination = (rank + 1) % size;
	sender = (rank - 1 + size) % size;

	this->x.resize(count_for_process);
	this->b.resize(count_for_process);
	for (int i = ibeg, j = 0; i < iend; i++, j++)
	{
		this->x[j] = x[i];
		this->b[j] = b[i];
	}
	norm_denominator = find_norm_b();
};

int Solving_Linear_Equations_parallel_second::find_norm_b()
{
	double sum_tmp = 0;
	double sum = 0;
	std::vector<double> tmp(count_for_process);
	for (int i = ibeg; i < iend; i++)
	{
		sum_tmp += b[i] * b[i];
	}
	MPI_Allreduce(&sum_tmp, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return sum;
}

double Solving_Linear_Equations_parallel_second::multiply_row_by_column(const std::vector<double>& row, const std::vector<double>& column, int offset, int count) const
{
	double result = 0.0;
	for (int i = offset, j = 0; j < count; ++i, ++j) {
		result += row[i] * column[j];
	}
	return result;
}

void Solving_Linear_Equations_parallel_second::proximity_function()
{
	std::vector<double> x_process(count_for_process, 0.0);
	std::vector<double> tmp(count_for_process);
	for (int k = 0; k < size; k++)
	{
		int block = (rank + k) % size;

		for (int i = ibeg, j = 0; i < iend; i++, j++)
		{
			x_process[j] += multiply_row_by_column(A[i], x, block * count_for_process, count_for_process);
		}

		if (size == 1) continue;
		
		MPI_Sendrecv(&x[0], count_for_process, MPI_DOUBLE, destination, 0,
			&tmp[0], count_for_process, MPI_DOUBLE, sender, MPI_ANY_TAG,
			MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		x = tmp;
	}
	
	for (int i = 0; i < count_for_process; i++)
	{
		x_process[i] = x[i] - ti * (x_process[i] - b[i]);
	}
	x = x_process;
}

bool Solving_Linear_Equations_parallel_second::accuracy_check(double epsilon)
{ 
	std::vector<double> result(count_for_process, 0.0);
	std::vector<double> tmp(count_for_process);
	for (int k = 0; k < size; k++)
	{
		int block = (rank + k) % size;

		for (int i = ibeg, j = 0; i < iend; i++, j++)
		{
			result[j] += multiply_row_by_column(A[i], x, block * count_for_process, count_for_process);
		}

		if (size == 1) continue;

		MPI_Sendrecv(&x[0], count_for_process, MPI_DOUBLE, destination, 0,
			&tmp[0], count_for_process, MPI_DOUBLE, sender, MPI_ANY_TAG,
			MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		x = tmp;
	}

	double part_norm = 0;

	for (int i = 0; i < count_for_process; i++)
	{
		part_norm += (result[i] - b[i]) * (result[i] - b[i]);
	}

	double norm_numerator;

	MPI_Allreduce(&part_norm, &norm_numerator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	return norm_numerator / norm_denominator < epsilon * epsilon ? true : false;
}

std::vector<double> Solving_Linear_Equations_parallel_second::build_res_vec()
{
	std::vector<double> res_vec(N);
	std::vector<double> tmp(count_for_process);
	for (int k = 0, i = 0; k < size; k++)
	{
		for (int j = 0; j < count_for_process; j++, i++)
		{
			res_vec[i] = x[j];
		}

		if (size == 1) continue;

		MPI_Sendrecv(&x[0], count_for_process, MPI_DOUBLE, destination, 0,
			&tmp[0], count_for_process, MPI_DOUBLE, sender, MPI_ANY_TAG,
			MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		x = tmp;
	}
	return res_vec;
}

std::vector<double> Solving_Linear_Equations_parallel_second::execute(double epsilon)
{
	while (accuracy_check(epsilon) == false)
	{
		proximity_function();
	}
	return build_res_vec();
}
