#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "Matrix.hpp"

#pragma region

Matrix::CommaInput::CommaInput(Matrix& matrix) : matrix(matrix), index(0) {
}

Matrix::CommaInput::CommaInput(Matrix& matrix, double input) : CommaInput(matrix) {
	(*this) << input;
}

Matrix::CommaInput::CommaInput(Matrix& matrix, Input input) : CommaInput(matrix) {
	(*this) << input;
}

Matrix::CommaInput& Matrix::CommaInput::operator <<(double input) {
	if (index < matrix.Size()) {
		matrix(index++) = input;
	}
	return *this;
}

Matrix::CommaInput& Matrix::CommaInput::operator ,(double input) {
	return (*this) << input;
}

Matrix::CommaInput& Matrix::CommaInput::operator <<(Input input) {
	switch (input)
	{
	case Matrix::NEXT_ENTRY:
		index++;
		break;
	case Matrix::END_ROW:
		int n = matrix.Column();
		index += n - (index % n);
		break;
	}
	return *this;
}

Matrix::CommaInput& Matrix::CommaInput::operator ,(Input input) {
	return (*this) << input;
}

Matrix IdentityMatrix(int dimension) {
	Matrix m(dimension, dimension);
	for (int i = 0; i < dimension; i++) {
		m(i, i) = 1.0;
	}
	return m;
}

#pragma endregion

#pragma region Constructor

Matrix::Matrix(int row, int column, double initValue) : row(row), column(column) {
	if (Row() < 1 || Column() < 1) {
		throw "Exception : Invailed matrix size.\n";
		exit(EXIT_FAILURE);
	}

	this->values = new double[Size()];
	if (!Values()) {
		throw "Exception : Can not allocate memory.\n";
		exit(EXIT_FAILURE);
	}

	for (long i = 0; i < Size(); i++) {
		(*this)(i) = initValue;
	}
}

Matrix::Matrix(int row, int column, double* initValues, long size) : Matrix(row, column) {
	if (size > Size()) {
		size = Size();
	}
	for (long i = 0; i < size; i++) {
		(*this)(i) = initValues[i];
	}
}

Matrix::Matrix(int row, int column, double* initValues) : Matrix(row, column, initValues, Size(row, column)) {
}

Matrix::Matrix(const Matrix& matrix) : Matrix(matrix.Row(), matrix.Column(), matrix.Values(), matrix.Size()) {
}

Matrix::~Matrix() {
	if (Values()) {
		delete[] Values();
	}
	if (commaInput) {
		delete commaInput;
	}
}

#pragma endregion

#pragma region Function

bool Matrix::EqualsType(const Matrix& matrix) const {
	return Row() == matrix.Row() && Column() == matrix.Column();
}

bool Matrix::IsSquare() const {
	return Row() == Column();
}

double Matrix::Abs() const {
	if (!IsSquare()) {
		throw "Exception : Not a square matrixe.\n";
		exit(EXIT_FAILURE);
	}

	Matrix m = UpperTriangular();
	double f = 1.0;
	for (int i = 0; i < Row(); i++) {
		f *= m(i, i);
	}
	return f;
}

double Matrix::Trace() const {
	if (!IsSquare()) {
		throw "Exception : Not a square matrixe.\n";
		exit(EXIT_FAILURE);
	}

	double f = 0.0;
	for (int i = 0; i < Row(); i++) {
		f += (*this)(i, i);
	}
	return f;
}

Matrix Matrix::UpperTriangular() const {
	if (!IsSquare()) {
		throw "Exception : Not a square matrixe.\n";
		exit(EXIT_FAILURE);
	}

	Matrix m = Clone();
	const int N = Row();
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			if (abs(m(i, i)) < abs(m(j, i))) {
				m.SwapRow(i, j);
			}
		}
		double aii = m(i, i);
		if (aii == 0.0) {
			throw "Exception : Can not calculate.\n";
			exit(EXIT_FAILURE);
		}
		
		for (int j = 0; j < N; j++) {
			if (i >= j) {
				continue;
			}
			double aji = m(j, i);
			for (int k = 0; k < N; k++) {
				m(j, k) -= m(i, k) * aji / aii;
			}
		}
	}
	return m;
}

Matrix Matrix::LowerTriangular() const {
	if (!IsSquare()) {
		throw "Exception : Not a square matrixe.\n";
		exit(EXIT_FAILURE);
	}

	Matrix m = Clone();
	const int N = Row();
	for (int i = N - 1; i >= 0; i--) {
		for (int j = i - 1; j >= 0; j--) {
			if (abs(m(i, i)) < abs(m(j, i))) {
				m.SwapRow(i, j);
			}
		}
		double aii = m(i, i);
		if (aii == 0.0) {
			throw "Exception : Can not calculate.\n";
			exit(EXIT_FAILURE);
		}
		
		for (int j = 0; j < N; j++) {
			if (i >= j) {
				continue;
			}
			double aji = m(j, i);
			for (int k = 0; k < N; k++) {
				m(j, k) -= m(i, k) * aji / aii;
			}
		}
	}
	return m;
}

Matrix Matrix::Minor(int row, int column, int rowSize, int columnSize) const {
	Matrix m(rowSize, columnSize);
	for (int i = 0; i < rowSize; i++) {
		for (int j = 0; j < columnSize; j++) {
			m(i, j) = (*this)(i + row, j + column);
		}
	}
	return m;
}

Matrix Matrix::MinorRemove(int row, int column) const {
	Matrix m(Row() - 1, Column() - 1);
	for (int i = 0; i < Row() - 1; i++) {
		for (int j = 0; j < Column() - 1; j++) {
			m(i, j) = (*this)(i < row ? i : i+1, j < column ? j : j+1);
		}
	}
	return m;
}

double Matrix::Cofactor(int row, int column) const {
	return pow(-1.0, row + column) * MinorRemove(row, column).Abs();
}

Matrix Matrix::Cofactor() const {
	Matrix m(Row(), Column());
	for (int i = 0; i < Row(); i++) {
		for (int j = 0; j < Column(); j++) {
			m(i, j) = Cofactor(i, j);
		}
	}
	return m;
}

Matrix Matrix::Transposed() const {
	Matrix m(Column(), Row());
	for (int i = 0; i < Row(); i++) {
		for (int j = 0; j < Column(); j++) {
			m(j, i) = (*this)(i, j);
		}
	}
	return m;
}

Matrix Matrix::Inverse() const {
	if (!IsSquare()) {
		throw "Exception : Not a square matrixe.\n";
		exit(EXIT_FAILURE);
	}

	Matrix m = Clone();
	const int N = Row();
	Matrix e = IdentityMatrix(N);
	for (int i = 0; i < N; i++) {
		for (int j = i + 1; j < N; j++) {
			if (abs(m(i, i)) < abs(m(j, i))) {
				m.SwapRow(i, j);
			}
		}
		double aii = m(i, i);
		if (aii == 0.0) {
			throw "Exception : Can not calculate.\n";
			exit(EXIT_FAILURE);
		}

		for (int j = 0; j < N; j++) {
			m(i, j) /= aii;
			e(i, j) /= aii;
		}
		for (int j = 0; j < N; j++) {
			if (i == j) {
				continue;
			}
			double aji = m(j, i);
			for (int k = 0; k < N; k++) {
				m(j, k) -= m(i, k) * aji;
				e(j, k) -= m(i, k) * aji;
			}
		}
	}
	delete &e;
	return m;
}

void Matrix::SwapRow(int a, int b) {
	for (int i = 0; i < Column(); i++) {
		std::swap((*this)(a, i), (*this)(b, i));
	}
}

void Matrix::SwapColumn(int a, int b) {
	for (int i = 0; i < Row(); i++) {
		std::swap((*this)(i, a), (*this)(i, b));
	}
}

#pragma endregion

#pragma region Operator

Matrix::CommaInput& Matrix::operator <<(double f) {
	if (commaInput) {
		delete commaInput;
	}
	commaInput = new CommaInput(*this, f);
	return *commaInput;
}

Matrix::CommaInput& Matrix::operator <<(Input i) {
	if (commaInput) {
		delete commaInput;
	}
	commaInput = new CommaInput(*this, i);
	return *commaInput;
}

std::ostream& operator <<(std::ostream& os, const Matrix& m) {
	for (int i = 0; i < m.Row(); i++) {
		for (int j = 0; j < m.Column(); j++) {
			os << m(i, j) << ",";
		}
		os << std::endl;
	}
	return os;
}

bool operator ==(const Matrix& a, const Matrix& b) {
	if (!a.EqualsType(b)) {
		return false;
	}
	for (long i = 0; i < a.Size(); i++) {
		if (a(i) != b(i)) {
			return false;
		}
	}
	return true;
}

bool operator !=(const Matrix& a, const Matrix& b) {
	return !(a == b);
}

Matrix& operator +=(Matrix& a, const Matrix& b) {
	if (!a.EqualsType(b)) {
		throw "Exception : Invailed matrix size.\n";
		exit(EXIT_FAILURE);
	}

	for (long i = 0; i < a.Size(); i++) {
		a(i) += b(i);
	}
	return a;
}

Matrix& operator -=(Matrix& a, const Matrix& b) {
	if (!a.EqualsType(b)) {
		throw "Exception : Invailed matrix size.\n";
		exit(EXIT_FAILURE);
	}

	for (long i = 0; i < a.Size(); i++) {
		a(i) -= b(i);
	}
	return a;
}

Matrix& operator *=(Matrix& a, double f) {
	for (long i = 0; i < a.Size(); i++) {
		a(i) *= f;
	}
	return a;
}

Matrix& operator /=(Matrix& a, double f) {
	return a *= (1.0 / f);
}

Matrix& operator +(const Matrix& a, const Matrix& b) {
	if (!a.EqualsType(b)) {
		throw "Exception : Invailed matrix size.\n";
		exit(EXIT_FAILURE);
	}

	Matrix *re = new Matrix(a.Row(), a.Column());
	for (long i = 0; i < re->Size(); i++) {
		(*re)(i) = a(i) + b(i);
	}
	return *re;
}

Matrix& operator -(const Matrix& a, const Matrix& b) {
	if (!a.EqualsType(b)) {
		throw "Exception : Invailed matrix size.\n";
		exit(EXIT_FAILURE);
	}

	Matrix *re = new Matrix(a.Row(), a.Column());
	for (long i = 0; i < re->Size(); i++) {
		(*re)(i) = a(i) - b(i);
	}
	return *re;
}

Matrix& operator *(const Matrix& a, const Matrix& b) {
	if (a.Column() != b.Row()) {
		throw "Exception : Invailed matrix size.\n";
		exit(EXIT_FAILURE);
	}

	Matrix *re = new Matrix(a.Row(), b.Column());
	for (int i = 0; i < re->Row(); i++) {
		for (int j = 0; j < re->Column(); j++) {
			register double f = 0.0;
			for (int k = 0; k < a.Column(); k++) {
				f += a(i, k) * b(k, j);
			}
			(*re)(i, j) = f;
		}
	}
	return *re;
}

Matrix& operator *(const Matrix& a, double f) {
	Matrix *re = new Matrix(a.Row(), a.Column());
	for (long i = 0; i < re->Size(); i++) {
		(*re)(i) = a(i) * f;
	}
	return *re;
}

Matrix& operator *(double f, const Matrix& a) {
	return a * f;
}

Matrix& operator /(const Matrix& a, double f) {
	return a * (1.0 / f);
}

Matrix& operator ^(const Matrix& a, int n) {
	if (n < 1) {
		throw "Exception : Invailed numeric value.\n";
		exit(EXIT_FAILURE);
	}
	Matrix *re = new Matrix(a);
	while(--n) {
		*re = *re * a;
	}
	return *re;
}

#pragma endregion
