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

class Solving_Linear_Equations_virtual
{
protected:
	double ti = ti_plus;
	static const double ti_plus;
	static const double ti_minus;
	int N;
	Matrix A;
	std::vector<double> b;
	std::vector<double> x;
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

class Solving_Linear_Equations_parallel_first : public Solving_Linear_Equations_virtual
{
private:
	std::vector<double> result_part;
	std::vector<double> result;

	int size, rank, ibeg, iend, count_for_process;
	void proximity_function() override;
	bool accuracy_check(double epsilon) override;
public:
	Solving_Linear_Equations_parallel_first(Matrix A, std::vector<double> x, std::vector<double> b);

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
	for (int i = 0; i < row.size(); i++)
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

		Solving_Linear_Equations_parallel_first solver(A, x, b);
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

