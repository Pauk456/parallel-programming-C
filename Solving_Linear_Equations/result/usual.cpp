#include<mpi.h>
#include<vector>
#include<iostream>
#include<chrono>
#include<cmath>

class Matrix
{
private:
	std::vector< std::vector<double> > matrix;
	int N;
public:
	Matrix(int N);

	int size() {
		return N;
	}
	const std::vector<double>& operator[](int index) const;
	std::vector<double>& operator[](int index);
};

Matrix::Matrix(int size) : N(size)
{
	matrix.resize(N, std::vector<double>(N, 1.0));
	for (size_t i = 0; i < matrix.size(); i++)
	{
		for (size_t j = 0; j < matrix.size(); j++)
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

class Solving_Linear_Equations_virtual
{
protected:
	double ti = ti_plus;
	static const double ti_plus;
	static const double ti_minus;
	Matrix A;
	int N;
	std::vector<double> x;
	std::vector<double> b;
	double norm_denominator = 0;

	virtual void proximity_function() = 0;
	virtual bool accuracy_check(double epsilon) = 0;
	double multiply_row_by_column(const std::vector<double>& row, const std::vector<double>& column) const;
	double find_norm(const std::vector<double>& row) const;
public:
	Solving_Linear_Equations_virtual(Matrix A) : A(A), N(A.size()) {};
	Solving_Linear_Equations_virtual(Matrix A, std::vector<double> x, std::vector<double> b);
	virtual ~Solving_Linear_Equations_virtual() = default;

	virtual void print_result() = 0;
	virtual void execute(double epsilon);
};

class Solving_Linear_Equations_usual : public Solving_Linear_Equations_virtual
{
private:
	void proximity_function() override;
	bool accuracy_check(double epsilon) override;
public:
	Solving_Linear_Equations_usual(Matrix A, std::vector<double> x, std::vector<double> b)
		: Solving_Linear_Equations_virtual(A, x, b) {};

	void print_result() override;
};

const double Solving_Linear_Equations_virtual::ti_plus = 0.00001;
const double Solving_Linear_Equations_virtual::ti_minus = -0.00001;

Solving_Linear_Equations_virtual::Solving_Linear_Equations_virtual(Matrix A, std::vector<double> x, std::vector<double> b)
	: A(A), N(A.size()), x(x), b(b)
{
	norm_denominator = find_norm(b);
}

double Solving_Linear_Equations_virtual::multiply_row_by_column(const std::vector<double>& row, const std::vector<double>& column) const
{
	double result = 0.0;
	for (int i = 0; i < N; ++i) {
		result += row[i] * column[i];
	}
	return result;
}

double Solving_Linear_Equations_virtual::find_norm(const std::vector<double>& row) const
{
	double result = 0.0;
	for (size_t i = 0; i < row.size(); i++)
	{
		result += row[i] * row[i];
	}
	result = std::pow(result, 0.5);
	return result;
}

void Solving_Linear_Equations_virtual::execute(double epsilon)
{
	while (accuracy_check(epsilon) == false)
	{
		proximity_function();
	}
}

void Solving_Linear_Equations_usual::proximity_function()
{
	std::vector<double> new_x(N);
	for (int i = 0; i < N; i++)
	{
		new_x[i] = multiply_row_by_column(A[i], x);
	}

	for (int i = 0; i < N; i++)
	{
		new_x[i] = x[i] - ti * (new_x[i] - b[i]);
	}
	x = new_x;
}

bool Solving_Linear_Equations_usual::accuracy_check(double epsilon)
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

void Solving_Linear_Equations_usual::print_result()
{
	std::cout << "Count of elements X = " << x.size() << std::endl;
	for (size_t i = 0; i < x.size(); i++)
	{
		std::cout << x[i] << ' ';
	}
	std::cout << std::endl;
}


int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int N = 1000;

	Matrix A(N);
	std::vector<double> x(N, 0);
	std::vector<double> b(N, N + 1);
	try {
		auto start_time = std::chrono::high_resolution_clock::now();

		Solving_Linear_Equations_usual solver(A, x, b);
		solver.execute(0.0001);

		auto end_time = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
		if (rank == 0)
		{
			std::cout << "Time passed: " << duration.count() << " micsec" << std::endl;
			solver.print_result();
		}
		
	}
	catch (const std::exception& e)
	{
		std::cout << e.what() << std::endl;
	}

	MPI_Finalize();
}
