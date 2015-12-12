#include "Matrix.h"
#include"GaussElimination.h"
#include<iostream>

#define MAX_SPACE 10

Matrix::Matrix(int r, int c)
{
	this->m_rows = r;
	this->m_column = c;
	this->m_isReleased = false;
	this->m_data = new float*[this->m_rows];
	for (int i = 0; i < this->m_rows; i++)
	{
		this->m_data[i] = new float[this->m_column];
	}

	for (int i = 0; i < this->m_rows; i++)
	{
		for (int j = 0; j < this->m_column; j++)
		{
			this->m_data[i][j] = 0.0f;
		}
	}
}

Matrix::Matrix(const Matrix & that)
{
	this->m_rows = that.RowSize();
	this->m_column = that.ColumnSize();
	this->m_data = new float*[this->m_rows];
	for (int i = 0; i < this->m_rows; i++)
	{
		this->m_data[i] = new float[this->m_column];
	}
	for (int i = 0; i < this->m_rows; i++)
	{
		for (int j = 0; j < this->m_column; j++)
		{
			this->m_data[i][j] = that.GetDataAt(i,j);
		}
	}
}

Matrix & Matrix::operator=(const Matrix & that)
{
	this->m_rows = that.RowSize();
	this->m_column = that.ColumnSize();
	this->m_data = new float*[this->m_rows];
	for (int i = 0; i < this->m_rows; i++)
	{
		this->m_data[i] = new float[this->m_column];
	}
	for (int i = 0; i < this->m_rows; i++)
	{
		for (int j = 0; j < this->m_column; j++)
		{
			this->m_data[i][j] = that.GetDataAt(i, j);
		}
	}
	return *this;
}

Matrix & Matrix::operator=(Matrix && that)
{
	this->m_rows = that.RowSize();
	this->m_column = that.ColumnSize();
	this->m_data = new float*[this->m_rows];
	for (int i = 0; i < this->m_rows; i++)
	{
		this->m_data[i] = new float[this->m_column];
	}
	for (int i = 0; i < this->m_rows; i++)
	{
		for (int j = 0; j < this->m_column; j++)
		{
			this->m_data[i][j] = that.GetDataAt(i, j);
		}
	}
	that.Release();
	return *this;
}

Matrix::~Matrix()
{
	if (!this->m_isReleased)
	Release();
}

Matrix::Matrix(Matrix && that)
{
	this->m_rows = that.RowSize();
	this->m_column = that.ColumnSize();
	this->m_data = new float*[this->m_rows];
	for (int i = 0; i < this->m_rows; i++)
	{
		this->m_data[i] = new float[this->m_column];
	}
	for (int i = 0; i < this->m_rows; i++)
	{
		for (int j = 0; j < this->m_column; j++)
		{
			this->m_data[i][j] = that.GetDataAt(i, j);
		}
	}
	that.Release();
}

void Matrix::UpdateData(int i, int j, float data)
{
	assert(i < this->m_rows && j < this->m_column);
	this->m_data[i][j] = data;
}

std::string Matrix::MatrixDataToString()
{
	std::string printText = "";
	for (int i = 0; i < this->m_rows; i++)
	{
		for (int j = 0; j < this->m_column; j++)
		{
			printText += std::to_string(this->m_data[i][j]);
			int len = strlen(std::to_string(this->m_data[i][j]).c_str());
			for (int i = 0; i < MAX_SPACE;i++)
				printText += " ";
			
		}
		printText += "\n";
	}

	return printText;
}

Matrix Matrix::Multiply(const Matrix & m)
{
	assert(this->m_column == m.RowSize());
	Matrix ret(this->m_rows,m.ColumnSize());

	for (int i = 0; i < this->m_rows; i++)
	{
		for (int j = 0; j < m.ColumnSize(); j++)
		{
			float sum = 0;
			for (int k = 0; k < this->m_column; k++)
			{
				sum += this->m_data[i][k] * m[k][j];
			}
			ret.UpdateData(i, j, sum);
		}
	}
	return ret;
}

void Matrix::Release()
{
	for (int i = 0; i < this->m_rows; i++)
		delete[] this->m_data[i];

	delete[] this->m_data;
	this->m_isReleased = true;
}

bool Matrix::operator==(const Matrix & rhs)
{
	for (int i = 0; i < rhs.RowSize(); i++)
	{
		for (int j = 0; j < rhs.ColumnSize(); j++)
		{
			if (this->m_data[i][j] != rhs.GetDataAt(i, j)) return false;
		}
	}
	return true;
}

const bool Matrix::operator==(const Matrix & rhs) const
{
	for (int i = 0; i < rhs.RowSize(); i++)
	{
		for (int j = 0; j < rhs.ColumnSize(); j++)
		{
			if (this->m_data[i][j] != rhs.GetDataAt(i, j)) return false;
		}
	}
	return true;
}

Matrix & Matrix::operator+(const Matrix & rhs)
{
	for (int i = 0; i < this->m_rows; i++)
	{
		for (int j = 0; j < this->m_column; j++)
		{
			this->m_data[i][j] = this->m_data[i][j] + rhs.GetDataAt(i, j);
		}
	}
	return *this;
}


Matrix & Matrix::operator-(const Matrix & rhs)
{
	for (int i = 0; i < this->m_rows; i++)
	{
		for (int j = 0; j < this->m_column; j++)
		{
			this->m_data[i][j] = this->m_data[i][j] - rhs.GetDataAt(i, j);
		}
	}
	return *this;
}

Matrix::Row Matrix::operator[](int i)
{
	assert(i < this->m_rows);
	Row r(this->GetRow(i),this->m_column);
	return r;
}
 
const Matrix::Row Matrix::operator[](int i) const
{
	assert(i < this->m_rows);
	Row r(this->GetRow(i),this->m_column);
	return r;
}

void Matrix::MultiplyRowBy(int rowNum, float multiplyNumber)
{
	for (int i = 0; i < this->m_column; i++)
	{
		this->m_data[rowNum][i] *= multiplyNumber;
	}
}

void Matrix::SwapRows(int row, int swapBy)
{
	float * temp = new float[this->m_column];
	temp = this->m_data[swapBy];
	this->m_data[swapBy] = this->m_data[row];
	this->m_data[row] = temp;
}
void Matrix::SwapColumns(int col, int swapBy)
{
	assert(col < this->m_column && swapBy < this->m_column);
	if (col == swapBy) return;
	float * temp = new float[this->m_rows];
	for (int i = 0; i < this->m_rows; i++)
	{
		temp[i] = this->m_data[i][col];
	}

	for (int i = 0; i < this->m_column; i++)
	{
		this->m_data[i][col] = this->m_data[i][swapBy];
	}

	for (int i = 0; i < this->m_column; i++)
	{
		this->m_data[i][swapBy] = temp[i];
	}
}


void Matrix::MultiplyRowByAndAddToRow(int row, float multiplyBy, int addRow)
{
	float * temp = new float[this->m_column];

	for (int i = 0; i < this->m_column; i++)
		temp[i] = this->m_data[row][i];

	for (int i = 0; i < this->m_column; i++)
		temp[i] *= multiplyBy;

	for (int i = 0; i < this->m_column; i++)
		temp[i] += this->m_data[addRow][i];

	this->m_data[addRow] = temp;
}

Matrix Matrix::Pow(int i)
{
	assert(this->m_column == this->m_rows);
	Matrix m(this->m_rows,this->m_column);
	m = *this;
	for (int j = 0; j < i; j++)
		m = m.Multiply(*this);
	return m;
}

void Matrix::DivideEachElementInRowBy(int row, float divideBy)
{
	assert(row < this->m_rows);

	for (int i = 0; i < this->m_column; i++)
	{
		this->m_data[row][i] = (this->m_data[row][i] / divideBy);
	}
}

void Matrix::DivideEachElementInColumnBy(int col, float divideBy)
{
	assert(col < this->m_column);

	for (int i = 0; i < this->m_rows; i++)
	{
		this->m_data[i][col] = (this->m_data[i][col] / divideBy);
	}
}

float Matrix::FindMaxInColumn(int col)
{
	assert(col < this->m_column);
	float ret = 0;
	for (int i = 0; i < this->m_rows; i++)
	{
		if (std::abs(this->m_data[i][col]>ret)) ret = std::abs(this->m_data[i][col]);
	}
	return ret;
}

int Matrix::FindMaxIndiceInColumn(int col,int from)
{
	assert(col < this->m_column && from<this->m_rows);
	float ret = -999999;
	int j = 0;
	for (int i = from; i < this->m_rows; i++)
	{
		if (std::abs(this->m_data[i][col])>ret)
		{
			ret = std::abs(this->m_data[i][col]);
			j = i;
		}
	}
	return j;
}

Matrix Matrix::Inverse()
{
	std::cout << "\n\nM:\n" << this->MatrixDataToString();
	assert(!IsSingular());
	std::cout << "\n\nM:\n" << this->MatrixDataToString();
	Matrix ret(this->m_rows, this->m_column);
	std::vector<float> answers;
	for (int i = 0; i < this->m_rows; i++)
	{
		Matrix b(this->m_rows, 1);
		b.SetRow(i, 1);
		GaussElimination ge(*this, b);
		ge.SetScalling(true);
		ge.SetCompletePivoting(true);
		answers = ge.Calculate();
		
		for (int j = this->m_column-1; j >=0; j--)
		{
			ret.UpdateData(j, i, answers[j]);
		}
		std::cout << "\n\n" << ret.MatrixDataToString() << "\n\n";
	}

	std::cout << "\n\n** Checking is inversed?*****\n\n";
	Matrix I = ret.Multiply(*this);
	std::cout << I.MatrixDataToString();
	std::cout << "\n\n** Finished *****\n\n";
	return ret;
}

void Matrix::SetRow(int row,float data)
{
	assert(row < this->m_rows);
	for (int i = 0; i < this->m_column; i++)
	{
		this->m_data[row][i] = data;
	}
}

float Matrix::Determinate()
{
	float det = 1.0f;
	Matrix detMat(*this);
	GaussElimination ge(detMat);
	ge.CalulateUpper();
	Matrix u = ge.GetUpperMatrix();

	for (int i = 0; i < this->m_column; i++)
	{
		det *= u[i][i];
	}
	
	for (auto const & member : this->m_determinantOperation)
	{
		det *= member;
	}

	for (int i = 0; i < ge.GetSignChangeCounter(); i++)
	{
		det *= -1;
	}
	return det;
}

void Matrix::ScaleMatrix()
{
	for (int i = 0; i < this->m_rows; i++)
	{
		float max = this->FindMaxInRow(i);
		if (max == 0) max = 1;
		this->DivideEachElementInRowBy(i, max);
		this->m_determinantOperation.push_back(max);
	}
}

float Matrix::FindMaxInRow(int row)
{
	int max = 0;
	assert(row < this->m_rows);
	float ret = -999999;
	for (int i = 0; i < this->m_column; i++)
	{
		if (std::abs(this->m_data[row][i])>ret)
		{
			ret = std::abs(this->m_data[row][i]);
		}
	}
	return ret;
}

bool Matrix::IsSingular()
{
	if (this->Determinate() == 0) return true;
	return false;
}

void Matrix::FindMaxInBlock(int i, int *reti, int *retj)
{
	float max = 0;
	for (int j = i; j < this->m_rows; j++)
	{
		for (int k = i; k < this->m_column; k++)
		{
			if (std::abs(this->m_data[j][k])>max)
			{
				max = std::abs(this->m_data[j][k]);
				*reti = j;
				*retj = k;
			}
		}
	}
}

void Matrix::LoadIdentity()
{
	for (int i = 0; i < this->m_rows; i++)
	{
		for (int j = 0; j < this->m_column; j++)
		{
			this->m_data[i][i] = 1;
		}
	}
}