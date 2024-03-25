#include<omp.h>
#include<vector>
#include<iostream>
#include<chrono>
#include <math.h>


const double t_plus = 0.0001;
const double t_minus = -0.0001;
const double epsilon = 0.0001;

class Matrix
{
private:
	std::vector< std::vector<double> > matrix;
	int N;

public:
	Matrix(int N)
	{

		N = N;
		matrix.resize(N, std::vector<double>(N, 1.0));
#pragma omp parallel for
		for (int i = 0; i < matrix.size(); i++)
		{
#pragma omp parallel for
			for (int j = 0; j < matrix.size(); j++)
			{
				if (i == j)
				{
					matrix[i][j] = 2.0;
				}
			}
		}
	}
	int size() {
		return N;
	}
	const std::vector<double>& operator[](int index) const;
	std::vector<double>& operator[](int index);
};

std::vector<double>& Matrix::operator[](int index) {
	return matrix[index];
}

const std::vector<double>& Matrix::operator[](int index) const {
	return matrix[index];
}


double find_norm(const std::vector<double>& row)
{
	double result = 0.0;
	for (int i = 0; i < row.size(); i++)
	{
		result += row[i] * row[i];
	}

	return result;
}

double multiply_row_by_column(const std::vector<double>& row, const std::vector<double>& column)
{

	double result = 0.0;
	for (int i = 0; i < row.size(); ++i) {
		result += row[i] * column[i];
	}
	return result;
}




int main(int argc, char* argv[])
{

	int N = 120;
	Matrix A(N);
	std::vector<double> b(N, N + 1);

	double b_norm = find_norm(b);

	std::vector<double> result(N);
	std::vector<double> x1(N, 2.0);
	std::vector<double> new_x(N);

	double ep = epsilon * epsilon;

	double res;

	double time_spent = 0.0;
	clock_t begin = clock();

	int size = omp_get_max_threads();
	if (N % size != 0)
	{
		printf("sadfsadfsfdaasdfasdf\n");
		return 1;
	}

	do {
		res = 0.0;


#pragma omp parallel for
		for (int i = 0; i < N; i++) {
			new_x[i] = multiply_row_by_column(A[i], x1);
		}
#pragma omp parallel for
		for (int i = 0; i < N; i++) {
			new_x[i] = x1[i] - t_plus * (new_x[i] - b[i]);
		}

		x1 = new_x;
#pragma omp parallel for
		for (int i = 0; i < N; i++) {
			result[i] = multiply_row_by_column(A[i], x1) - b[i];
		}
#pragma omp parallel for reduction(+:res)
		for (int i = 0; i < N; i++) {
			res += result[i] * result[i];
		}

#pragma omp for reduction(+:b_norm)
		for (int i = 0; i < b.size(); i++) {
			b_norm += b[i] * b[i];
		}


	} while (!(res / b_norm < ep));

	clock_t end = clock();
	time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
	printf("time spend: %f \n", time_spent);

	return 0;
}
