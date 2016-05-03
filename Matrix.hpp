#pragma once

#pragma region Extern

class Matrix
{
public:
	enum Input {
		NEXT_ENTRY,
		END_ROW
	};
private:
	class CommaInput {
	public:
		CommaInput(Matrix& matrix, double input);
		CommaInput(Matrix& matrix, Input input);
		CommaInput& operator <<(double);
		CommaInput& operator ,(double);
		CommaInput& operator <<(Input);
		CommaInput& operator ,(Input);
	private:
		CommaInput(Matrix& matrix);
		long index;
		Matrix& matrix;
	};
public:
	Matrix(int row, int column, double initValue = 0.0);
	Matrix(int row, int column, double* initValues, long size);
	Matrix(int row, int column, double* initValues);
	Matrix(const Matrix& matrix);
	~Matrix();

	inline Matrix Clone() const;
	inline int Row() const;
	inline int Column() const;
	inline long Size() const;
	inline long Size(int row, int column) const;

	bool EqualsType(const Matrix& matrix) const;
	bool IsSquare() const;
	double Abs() const;
	double Trace() const;

	Matrix UpperTriangular() const;
	Matrix LowerTriangular() const;
	Matrix Minor(int row, int column, int rowSize, int columnSize) const;
	Matrix MinorRemove(int row, int column) const;
	double Cofactor(int row, int column) const;
	Matrix Cofactor() const;
	Matrix Transposed() const;
	Matrix Inverse() const;

	inline double& operator ()(long) const;
	inline double& operator ()(int, int) const;

	CommaInput& operator <<(double);
	CommaInput& operator <<(Input);
	friend std::ostream& operator <<(std::ostream&, const Matrix&);

	inline Matrix& operator +();
	inline Matrix& operator -();

	friend bool operator ==(const Matrix&, const Matrix&);
	friend bool operator !=(const Matrix&, const Matrix&);
	friend Matrix& operator +=(Matrix&, const Matrix&);
	friend Matrix& operator -=(Matrix&, const Matrix&);
	friend Matrix& operator *=(Matrix&, double);
	friend Matrix& operator /=(Matrix&, double);
	friend Matrix& operator +(const Matrix&, const Matrix&);
	friend Matrix& operator -(const Matrix&, const Matrix&);
	friend Matrix& operator *(const Matrix&, const Matrix&);
	friend Matrix& operator *(const Matrix&, double);
	friend Matrix& operator *(double, const Matrix&);
	friend Matrix& operator /(const Matrix&, double);
	friend Matrix& operator ^(const Matrix&, int);
protected:
	inline double* Values() const;

	inline void Row(int row);
	inline void Column(int column);

	void SwapRow(int a, int b);
	void SwapColumn(int a, int b);
private:
	int row;
	int column;
	double *values;

	CommaInput* commaInput;
};

Matrix IdentityMatrix(int dimension);

#pragma endregion

#pragma region Getter

Matrix Matrix::Clone() const {
	return *(new Matrix(*this));
}

int Matrix::Row() const {
	return this->row;
}

int Matrix::Column() const {
	return this->column;
}

long Matrix::Size() const {
	return Size(Row(), Column());
}

long Matrix::Size(int row, int column) const {
	return (long)row * column;
}

double* Matrix::Values() const {
	return this->values;
}

#pragma endregion

#pragma region Setter

void Matrix::Row(int row) {
	this->row = row;
}

void Matrix::Column(int column) {
	this->column = column;
}

#pragma endregion

#pragma region Operator

double& Matrix::operator ()(long index) const {
	return this->values[index];
}

double& Matrix::operator ()(int row, int column) const {
	return (*this)((long)Column() * row + column);
}

Matrix& Matrix::operator +() {
	return Clone();
}

Matrix& Matrix::operator -() {
	return (*this) * -1.0;
}

#pragma endregion
