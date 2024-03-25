#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

typedef struct {
	int rows;
	int cols;
	double** data;
} Matrix;

int Matrix_create(Matrix *Matrix)
{
	Matrix->data = malloc(sizeof(double*) * Matrix->rows);
	if (!Matrix->data)
	{
		return ENOMEM;
	}
	for (int i = 0; i < Matrix->rows; i++)
	{
		Matrix->data[i] = malloc(sizeof(double) * Matrix->cols);
		if (!Matrix->data[i])
		{
			for (int j = 0; j < i; j++)
			{
				free(Matrix->data[j]);
			}
			free(Matrix->data);
			return ENOMEM;
		}
	}
	return 0;
}

void Matrix_free(Matrix *Matrix)
{
	for (int i = 0; i < Matrix->rows; i++)
	{
		free(Matrix->data[i]);
	}
	free(Matrix->data);
}

void Matrix_initialization(Matrix* Matrix)
{
	for (int i = 0; i < Matrix->rows; i++) {
		for (int j = 0; j < Matrix->cols; j++) {
			if (i == j) {
				Matrix->data[i][j] = 2.0;
			}
			else {
				Matrix->data[i][j] = 1.0;
			}
		}
	}
}

void forward_step1(MPI_Comm commGrib)
{
	
}

void forward_everyone(MPI_Comm commGrib)
{
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	int N1, N2, N3, N4;

	N1 = 100;
	N2 = N3 = 100;
	N4 = 100;
	
	int dims[2] = { 0,0 }, coords[2], reorder = 1;
	int size, rank, sizey, sizex, ranky, rankx;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Dims_create(size, 2, dims);
	sizey = dims[0]; sizex = dims[1];
	MPI_Comm commGrib;
	int periods[2] = { 0,0 }, reorder = 1;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &commGrib);

	MPI_Datatype row_type, col_type;
	MPI_Type_contiguous(N1, MPI_DOUBLE, &row_type);
	MPI_Type_vector(N, 1, N, MPI_DOUBLE, &col_type);
	MPI_Type_commit(&row_type);
	MPI_Type_commit(&col_type);


	Matrix Matrix1;
	Matrix1.rows = N1;
	Matrix1.cols = N2;
	if (rank == 0)
	{
		Matrix_create(&Matrix1);
		Matrix_initialization(&Matrix1);
	}
	


	if (rank == 0)
	{
		Matrix_free(&Matrix1);
	}
	MPI_Finalize();
}