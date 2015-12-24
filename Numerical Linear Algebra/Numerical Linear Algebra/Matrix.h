#pragma once
#include<string>
#include<assert.h>
#include<vector>

class GaussElimination;

class Matrix
{
public:

	class Row
	{
	public:
		Row(float * data,int size)
		{
			this->m_size = size;
			this->m_data = data;
		}
		Row(const float * data,int size)
		{
			this->m_size = size;
			this->m_data = new float[size];
			for (int i = 0; i < size; i++)
			{
				this->m_data[i] = data[i];
			}
			
		}

		float &operator [](const int j)
		{
			assert(j < this->m_size);
			return this->m_data[j];
		}


		const float &operator [](const int j) const
		{
			assert(j < this->m_size);
			return this->m_data[j];
		}
	private:
		float * m_data;
		int m_size;
	
	};

public:

	Matrix(int r,int c);
	Matrix() {}
	Matrix(const Matrix & that);
	Matrix(Matrix && that);
	~Matrix();
	Matrix & operator =(const Matrix & that);
	Matrix & operator =(Matrix &&that);

public:
	void UpdateData(int i, int j, float data);
	std::string MatrixDataToString();
	Matrix Multiply(const Matrix & m);
	void Release();
	void MultiplyRowBy(int rowNum, float multiplyNumber);
	void SwapRows(int row, int swapBy);
	void SwapColumns(int col, int swapBy);
	void MultiplyRowByAndAddToRow(int row, float multiplyBy, int addRow);
	void DivideEachElementInRowBy(int row, float divideBy);
	void DivideEachElementInColumnBy(int col, float divideBy);
	Matrix Pow(int i);
	float FindMaxInColumn(int col);
	float FindMaxInRow(int row);
	int FindMaxIndiceInColumn(int col,int from);
	Matrix Inverse();
	void SetRow(int row, float data);
	float Determinate();
	void ScaleMatrix();
	bool IsSingular();
	void FindMaxInBlock(int i, int *reti, int *retj);
	void LoadIdentity();
	Matrix Transpose() const;
	bool IsIdentity();
	bool IsOrthogonal();
	float MultiplyRowByColumnOfMatrixA(int row, const Matrix &A, int col);
	void QRFactorization();
	void LoadIdentityUpperTriangle();
	void FillMatrix(float fillWith);
	void Divide(float divideBy);
	float SumOfRowAfterIndicePow2(int row,int col);
	Matrix GetColumnInMatrix(int col);
	void SetColumn(int col,const Matrix &mat);
	void SetSubMatrix(const Matrix &mat, int fromRow, int fromCol);
	void SetColumnOfMatrixFromRow(const Matrix &mat, int col,int row);
	float  GetEigenvalues(float precision);
	void GetCharacteristicPolynomial();
	float Trace();
	float GetInfiniteNorm();

public:
	int ColumnSize() { return this->m_column; }
	const int ColumnSize() const { return this->m_column; }
	int RowSize() { return this->m_rows; }
	const int RowSize() const { return this->m_rows; }
	const float GetDataAt(int i, int j) const { return this->m_data[i][j]; }
	float GetDataAt(int i, int j)  { return this->m_data[i][j]; }
	const float * GetRow(int i) const { return this->m_data[i]; }
	float * GetRow(int i) { return this->m_data[i]; }

public:
	bool operator ==(const Matrix & that);
	const bool operator ==(const Matrix & that) const;

	Matrix & operator +(const Matrix & rhs);
	Matrix & operator -(const Matrix & rhs);

	void operator /(const float divideBy);
	const void operator /(const float divideBy) const;

	Matrix& operator *(const float i);
	const Matrix& operator *(const float i) const;

	Matrix operator *(const Matrix & rhs);

	Matrix& operator *=(const float i);
	const Matrix& operator *=(const float i) const;

	Matrix& operator /=(const float divideBy);
	const Matrix& operator /=(const float divideBy) const;

	Row operator[](int i);
	const Row operator [](int i) const;

	
private:
	int m_column, m_rows;
	bool m_isReleased;
	float ** m_data;
	std::vector<float> m_determinantOperation;

};
