#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

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

static double multiply_row_by_column(double* row, double* column, int N)
{
	double result = 0.0;
	for (int i = 0; i < N; ++i) {
		result += row[i] * column[i];
	}
	return result;
}

static double find_norm(double* vector, int N,
	int ibeg, int iend, double general)
{
#pragma omp barrier
#pragma omp single
	{
		general = 0;
	}
	

	double partResult = 0.0;
	for (size_t i = ibeg; i < iend; i++)
	{
		partResult += vector[i] * vector[i];
	}
#pragma omp critical
	{
		general += partResult;
	}
#pragma omp barrier
#pragma omp single
	{
		general = sqrt(general);
	}
	return general;
}

static int Solving_Linear_Equations_First(double** A, double* x, double* b, double epsilon, int N)
{
	double* tempArray = malloc(sizeof(double) * N);
	if (!tempArray)
	{
		return ENOMEM;
	}
	double norm_denominator;
	double norm_numerator;
	double general = 0.0;

	int size = omp_get_max_threads();
	if (N % size != 0)
	{
		printf("sssssssssssssssssss\n");
		exit(1);
	}
	int count_for_process = N / size;

	do
	{
	#pragma omp parallel
		{

			norm_denominator = find_norm(b, N, ibeg, iend, general);
	
			// accuracy cheack
			for (int i = ibeg; i < iend; i++)
			{
				tempArray[i] = multiply_row_by_column(A[i], x, N);
			}

			for (int i = ibeg; i < iend; i++)
			{
				tempArray[i] = (tempArray[i] - b[i]);
			}
			norm_numerator = find_norm(b, N, ibeg, iend, general);
			// accuracy cheack

			// proximity
			for (int i = ibeg; i < iend; i++)
			{
				tempArray[i] = multiply_row_by_column(A[i], x, N);
			}

			for (int i = ibeg; i < iend; i++)
			{
				x[i] -= ti * (tempArray[i] - b[i]);
			}
			// proximity
	#pragma omp barrier
		}
	} while (!(norm_numerator / norm_denominator < epsilon));

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
	int N = 12;
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