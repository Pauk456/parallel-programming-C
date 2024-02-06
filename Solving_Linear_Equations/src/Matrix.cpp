#include "Matrix.h"

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
