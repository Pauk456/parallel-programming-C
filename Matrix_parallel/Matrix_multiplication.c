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

void grib_initialization(MPI_Comm commGrib, MPI_Datatype row_type, MPI_Datatype col_type, Matrix matrix1, Matrix matrix2)
{
	int dims[2] = { 0,0 }, coords[2], reorder = 1;
	int size, rank, sizey, sizex, ranky, rankx;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Dims_create(size, 2, dims);
	sizey = dims[0]; sizex = dims[1];
	
	int periods[2] = { 0,0 }, reorder = 1;

	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &commGrib);
	
	MPI_Type_contiguous(matrix1.rows, MPI_DOUBLE, &row_type);
	MPI_Type_vector(matrix2.cols, 1, matrix2.rows, MPI_DOUBLE, &col_type);
	MPI_Type_commit(&row_type);
	MPI_Type_commit(&col_type);
}

void colm_row_initialization(MPI_Comm commGrib, MPI_Comm commCol, MPI_Comm commRow)
{
	int remains_dims[2];
	remains_dims[0] = true;
	remains_dims[1] = false;
	MPI_Cart_sub(commGrib, remains_dims, commColm); //  Создание столбцов
	MPI_Scatterv(); // Рассылает?
}

void forward_step1(MPI_Comm commGrib, int size, int rank)
{

}

void forward_everyone(MPI_Comm commGrib, MPI_Datatype row_type, MPI_Datatype col_type)
{
	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	forward_step1(commGrib, rank, size);
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	int N1, N2, N3, N4;

	N1 = 100;
	N2 = N3 = 100;
	N4 = 100;
	
	Matrix Matrix1;
	Matrix1.rows = N1;
	Matrix1.cols = N2;
	if (rank == 0)
	{
		Matrix_create(&Matrix1);
		Matrix_initialization(&Matrix1);
	}

	MPI_Comm commGrib;
	MPI_Datatype row_type, col_type;
	grib_initialization(commGrib, rowType, colType, matrix1, matrix2);
	forward_everyone(commGrib, row_type, col_type);
	


	if (rank == 0)
	{
		Matrix_free(&Matrix1);
	}
	MPI_Finalize();
}