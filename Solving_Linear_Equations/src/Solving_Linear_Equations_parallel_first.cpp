#include "Solving_Linear_Equations.h"

#define FIRST_THREAD 0

std::chrono::microseconds durationGLOBAL = std::chrono::microseconds(0);

// Возможно можно увеличить производительность если process_collector = size - 1 
void Solving_Linear_Equations_parallel_first::proximity_function()
{ 
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int count_for_process = ceil((double)N / size);
	int ibeg = count_for_process * rank;
	int iend = count_for_process * (rank + 1) > N ? N : count_for_process * (rank + 1);

	std::vector<double> x_process(iend - ibeg);
	for (int i = ibeg, j = 0; i < iend; i++, j++) // ЫЫЫЫЫ
	{
		x_process[j] = multiply_row_by_column(A[i], x);
	}
	for (int i = ibeg, j = 0; i < iend; i++, j++)
	{
		x_process[j] = x[i] - ti * (x_process[j] - b[i]);
	}

	int process_collector = 0;
	if (rank == FIRST_THREAD)
	{
		auto start_time = std::chrono::high_resolution_clock::now();
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
		auto end_time = std::chrono::high_resolution_clock::now();
		durationGLOBAL += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
	}
	else
	{
		MPI_Ssend(&x_process[0], x_process.size(), MPI_DOUBLE, process_collector, 0, MPI_COMM_WORLD);
		MPI_Recv(&x[0], x.size(), MPI_DOUBLE, process_collector, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
}

bool Solving_Linear_Equations_parallel_first::accuracy_check(double epsilon) const
{ // Можно ли сделать так чтобы первый процесс который до сюда дошел стал process_collectorом?
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int count_for_process = ceil((double)N / size);
	int ibeg = count_for_process * rank;
	int iend = count_for_process * (rank + 1) > N ? N : count_for_process * (rank + 1);

	std::vector<double> result(N);
	for (int i = ibeg; i < iend; i++)
	{
		result[i] = multiply_row_by_column(A[i], x);
		result[i] = (result[i] - b[i]);
	}
	
	int process_collector = 0;
	if (rank == process_collector)
	{
		for (int i = 1; i < size; i++)
		{
			MPI_Recv(&result[i * count_for_process], count_for_process, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		for (int i = 1; i < size; i++)
		{
			MPI_Ssend(&result[0], result.size(), MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
	}
	else
	{
		MPI_Ssend(&result[ibeg], iend - ibeg, MPI_DOUBLE, process_collector, 0, MPI_COMM_WORLD);
		MPI_Recv(&result[0], result.size(), MPI_DOUBLE, process_collector, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	double norm_numerator = find_norm(result);
	return norm_numerator / norm_denominator < epsilon* epsilon ? true : false;
}
