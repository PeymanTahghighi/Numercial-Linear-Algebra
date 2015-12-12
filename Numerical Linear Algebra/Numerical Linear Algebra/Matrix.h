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

		float operator [](int j)
		{
			assert(j < this->m_size);
			return this->m_data[j];
		}

		const float operator [](int j) const
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

	Row operator[](int i);
	const Row operator [](int i) const;
	
private:
	int m_column, m_rows;
	bool m_isReleased;
	float ** m_data;
	std::vector<float> m_determinantOperation;

};
