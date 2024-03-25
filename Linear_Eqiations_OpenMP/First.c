#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

static double ti = 0.00001;

static void initialization_A(double** A, int N)
{
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			if (i == j)
			{
				A[i][j] = 2.0;
			}
			else
			{
				A[i][j] = 1;
			}
		}
	}
}

static void initialization_X(double* x, int N)
{
	for (size_t i = 0; i < N; i++)
	{
		x[i] = 0;
	}
}

static void initialization_B(double* b, int N)
{
	for (size_t i = 0; i < N; i++)
	{
		b[i] = N + 1;
	}
}

static double find_norm(double* vector, int N)
{
	double result = 0.0;
	#pragma omp parallel for
	for (size_t i = 0; i < N; i++)
	{
		result += vector[i] * vector[i];
	}
	result = sqrt(result);
	return result;
}

static double multiply_row_by_column(double* row, double* column, int N)
{
	double result = 0.0;
	#pragma omp parallel for
	for (int i = 0; i < N; ++i) {
		result += row[i] * column[i];
	}
	return result;
}

static int accuracy_check(double** A, double* x, double* b, double epsilon, double* tempArray, double norm_denominator, int N)
{
	#pragma omp parallel for
	for (int i = 0; i < N; i++)
	{
		tempArray[i] = multiply_row_by_column(A[i], x, N);
	}

	#pragma omp parallel for
	for (int i = 0; i < N; i++)
	{
		tempArray[i] = (tempArray[i] - b[i]);
	}

	double norm_numerator = find_norm(tempArray, N);
	return norm_numerator / norm_denominator < epsilon;
}

static void proximity_function(double** A, double* x, double* b, double* tempArray, int N)
{
	#pragma omp parallel for
	for (int i = 0; i < N; i++)
	{
		tempArray[i] = multiply_row_by_column(A[i], x, N);
	}

	#pragma omp parallel for
	for (int i = 0; i < N; i++)
	{
		x[i] -= ti * (tempArray[i] - b[i]);
	}
}

static int Solving_Linear_Equations_First(double** A, double* x, double* b, double epsilon, int N)
{
	double* tempArray = malloc(sizeof(double) * N);
	if (!tempArray)
	{
		return ENOMEM;
	}

	double norm_denominator = find_norm(b, N);

	while (!accuracy_check(A, x, b, epsilon, tempArray, norm_denominator, N))
	{
		proximity_function(A, x, b, tempArray, N);
	}
	free(tempArray);
}

static void print_vector(double* x, int N)
{
	printf("Count of elements X = %d\n", N);
	for (size_t i = 0; i < N; i++)
	{
		printf("%lf ", x[i]);
	}
	printf("\n");
}

int main()
{
	int N = 100;
	double epsilon = 0.0001;

	double* x = malloc(sizeof(double) * N);

	if (!x)
	{
		printf("error1\n");
		return ENOMEM;
	}

	double* b = malloc(sizeof(double) * N);

	if (!b)
	{
		free(x);
		printf("error2\n");
		return ENOMEM;
	}

	double** A = malloc(N * sizeof(double*));
	if (!A)
	{
		free(x);
		free(b);
		printf("error3\n");
		return ENOMEM;
	}

	for (size_t i = 0; i < N; i++)
	{
		A[i] = malloc(N * sizeof(double));
		if (!A[i])
		{
			for (size_t j = 0; j < i; j++)
			{
				free(A[j]);
			}
			free(A);
			free(x);
			free(b);
			printf("error4\n");
			return ENOMEM;
		}
	}

	initialization_X(x, N);
	initialization_B(b, N);
	initialization_A(A, N);

	double time_spent = 0.0;
	clock_t begin = clock();

	if (Solving_Linear_Equations_First(A, x, b, epsilon, N) == ENOMEM)
	{
		free(x);
		free(b);
		free(A);
		printf("error5\n");
		return ENOMEM;
	}

	clock_t end = clock();
	time_spent += (double)(end - begin) / CLOCKS_PER_SEC;
	printf("time spend: %f \n", time_spent);

	print_vector(x, N);

	free(x);
	free(b);
	for (int i = 0; i < N; i++)
	{
		free(A[i]);
	}
	free(A);
	return 0;
}