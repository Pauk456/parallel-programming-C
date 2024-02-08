#pragma once
#include<vector>

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
