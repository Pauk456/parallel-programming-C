#include "Solving_Linear_Equations.h"


const double Solving_Linear_Equations::ti_plus = 0.01;
const double Solving_Linear_Equations::ti_minus = -0.01;

Matrix::Matrix(int size) : N(size)
{
	matrix.resize(N, std::vector<double>(N, 1.0));

	for (int i = 0; i < matrix.size(); i++)
	{
		for (int j = 0; j < matrix.size(); j++)
		{
			if (i == j)
			{
				matrix[i][j] = 2.0;
			}
		}
	}
}

std::vector<double>& Matrix::operator[](int index) {
	return matrix[index];
}

const std::vector<double>& Matrix::operator[](int index) const {
	return matrix[index];
}

Solving_Linear_Equations::Solving_Linear_Equations(Matrix A, std::vector<double> x, std::vector<double> b) 
	: A(A), N(A.size()), x(x), b(b)
{
	norm_denominator = find_norm(b);
}

double Solving_Linear_Equations::multiply_row_by_column(const std::vector<double>& row, const std::vector<double>& column) const{
	double result = 0.0;
	for (int i = 0; i < N; ++i) {
		result += row[i] * column[i];
	}
	return result;
}

double Solving_Linear_Equations::find_norm(const std::vector<double>& row) const
{
	double result = 0.0;
	for (int i = 0; i < row.size(); i++) 
	{
		result += row[i] * row[i];
	}
	result = std::pow(result, 0.5);
	return result;
}

void Solving_Linear_Equations::proximity_function()
{
	/*static const int process_root = 0;
	int process_size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &process_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == process_root)
	{
		
	}
	else
	{

		MPI_Gather(,1,MPI_DOUBLE, nullptr, 0, MPI_DOUBLE, process_root, MPI_COMM_WORLD);
	}*/
	std::vector<double> new_x(N);
	for (int i = 0; i < N; i++)
	{
		new_x[i] = multiply_row_by_column(A[i], x);
	}

	for (int i = 0; i < N; i++)
	{
		new_x[i] = x[i] - ti_plus * (new_x[i] - b[i]);
	}

	x = new_x;
}

bool Solving_Linear_Equations::accuracy_check(double epsilon) const
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

std::vector<double> Solving_Linear_Equations::execute(double epsilon)
{
	while (accuracy_check(epsilon) == false)
	{
		proximity_function();
	}
	return x;
}
