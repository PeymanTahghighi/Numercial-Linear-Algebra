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
	
	{
		this->m_data = new float*[this->m_rows];
		for (int i = 0; i < this->m_rows; i++)
		{
			this->m_data[i] = new float[this->m_column];
		}
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



Matrix Matrix::Transpose() const
{
	Matrix ret(this->m_column, this->m_rows);
	for (int i = 0; i < this->m_rows; i++)
	{
		for (int j = 0; j < this->m_column; j++)
		{
			ret[j][i] = this->m_data[i][j];
		}
	}
	return ret;
}

bool Matrix::IsIdentity()
{
	if (this->m_column != this->m_rows) return false;
	for (int i = 0; i < this->m_rows; i++)
	{
		if (this->m_data[i][i] != 1) return false;
	}
	return true;
}

bool Matrix::IsOrthogonal()
{
	Matrix mt(*this);
	Matrix m(*this);
	mt.Transpose();
	if (mt.Multiply(m).IsIdentity()) return true;
	return false;
}

float Matrix::MultiplyRowByColumnOfMatrixA(int row, const Matrix &A, int col)
{
	assert(this->m_rows == A.ColumnSize());
	float ret = 0;
	for (int i = 0; i < this->m_rows; i++)
	{
		ret += this->m_data[row][i] * A[i][col];
	}
	return ret;
}

void Matrix::QRFactorization()
{
	int k = 0;
	Matrix R(*this);
	std::vector<Matrix> Qs;
	while (k < this->m_column-1)
	{
		int sign = 0;
		R[k][k] > 0 ? sign = 1 : sign = -1;
		Matrix H(this->m_rows, this->m_rows);
		H.LoadIdentity();
		Matrix tmpH(this->m_rows - k, this->m_rows - k);
		Matrix I(this->m_rows - k, this->m_rows - k);
		Matrix vk(this->m_rows - k, 1);
		Matrix ek(this->m_rows - k, 1);
		Matrix vkT(1, this->m_rows - k);
		Matrix xk(this->m_rows - k, 1);
		H.LoadIdentity();
		
		xk.SetColumnOfMatrixFromRow(*this, k, k);
		//std::cout << "*****\n\nxk:\n" << xk.MatrixDataToString();
		float xknorm = sqrtf(xk.SumOfRowAfterIndicePow2(0, 0));
		
		ek[0][0] = 1.0f;
		vk = xk - ek*(xknorm*sign);
		//std::cout << "*****\n\vk:\n" << vk.MatrixDataToString();
		
		vkT = vk.Transpose();
		//std::cout << "*****\n\vkt:\n" << vkT.MatrixDataToString();
		Matrix vktV = vkT.Multiply(vk);
		//std::cout << "*****\n\vktv:\n" << vktV.MatrixDataToString();
		float c = 2 / vktV[0][0];
		
		I.LoadIdentity();
		Matrix vvt(this->m_rows - k, this->m_rows - k);
		//std::cout << "*****\n\vk:\n" << vk.MatrixDataToString();
		//std::cout << "*****\n\vkt:\n" << vkT.MatrixDataToString();
		vvt = vk.Multiply(vkT);
		//std::cout << "*****\n\vvt:\n" << vvt.MatrixDataToString();
		vvt *= c;
		tmpH = I - vvt;
		//std::cout << "*****\n\TMPH:\n" << tmpH.MatrixDataToString();
		H.SetSubMatrix(tmpH, k, k);
		Qs.push_back(H);
		//std::cout << "*****\n\H:\n" << H.MatrixDataToString();
		*this = H.Multiply(*this);
		//std::cout << "*****\n\new Matrix\n" << this->MatrixDataToString();
		k++;
	}
	Matrix tmp;
	tmp = *this;
	*this = R;
	R=tmp;
	std::cout << "\n\n******\nMATRIX R:\n\n" << R.MatrixDataToString();

	Matrix Q(this->m_rows, this->m_rows);
	if (Qs.size() == 1)
	{
		Q = Qs[0];
	}

	for (unsigned int i = 0; i < Qs.size() - 1; i++)
	{

		Q = Qs[i]*Qs[i + 1];
	}
	std::cout << "\n\n******\nMATRIX Q:\n\n" << Q.MatrixDataToString();

	std::cout << "***/n/n Checking:" << std::endl << std::endl << std::endl << Q.Multiply(R).MatrixDataToString();

}

void Matrix::LoadIdentityUpperTriangle()
{
	int k = 1;
	for (int i = 0; i < this->m_rows; i++)
	{
		for (int j = 0; j < k; j++)
		{
			this->m_data[j][i] = 1;
		}
		k++;
	}
}

void Matrix::FillMatrix(float fillWith)
{
	for (int i = 0; i < this->m_rows; i++)
	{
		for (int j = 0; j < this->m_column; j++)
		{
			this->m_data[i][j] = fillWith;
		}
	}
}

void Matrix::Divide(float divideBy)
{
	for (int i = 0; i < this->m_rows; i++)
	{
		for (int j = 0; j < this->m_column; j++)
		{
			this->m_data[i][j] /= divideBy;
		}
	}
}

void Matrix::operator/(float divideBy)
{
	for (int i = 0; i < this->m_rows; i++)
	{
		for (int j = 0; j < this->m_column; j++)
		{
			this->m_data[i][j] /= divideBy;
		}
	}
}
const void Matrix::operator/(float divideBy) const
{
	for (int i = 0; i < this->m_rows; i++)
	{
		for (int j = 0; j < this->m_column; j++)
		{
			this->m_data[i][j] /= divideBy;
		}
	}
}

Matrix& Matrix::operator/=(float divideBy)
{
	for (int i = 0; i < this->m_rows; i++)
	{
		for (int j = 0; j < this->m_column; j++)
		{
			this->m_data[i][j] /= divideBy;
		}
	}
	return *this;
}
const Matrix& Matrix::operator/=(float divideBy) const
{
	for (int i = 0; i < this->m_rows; i++)
	{
		for (int j = 0; j < this->m_column; j++)
		{
			this->m_data[i][j] /= divideBy;
		}
	}
	return *this;
}

float Matrix::SumOfRowAfterIndicePow2(int row, int col)
{
	assert(row <= this->m_rows && col<=this->m_column);
	float ret = 0;
	for (int i = row; i < this->m_rows; i++)
	{
		ret += pow(this->m_data[i][col],2);
	}
	return ret;
}

Matrix Matrix::GetColumnInMatrix(int col)
{
	assert(col < this->m_column);
	Matrix ret(this->m_rows, 1);
	for (int i = 0; i < this->m_rows; i++)
	{
		ret[i][0] = this->m_data[i][col];
	}
	return ret;
}

void Matrix::SetColumn(int col,const Matrix & mat)
{
	assert(mat.ColumnSize() == 1 && col<this->m_column);

	for (int i = 0; i < this->m_rows; i++)
	{
		this->m_data[i][col] = mat[i][0];
	}
}

Matrix& Matrix::operator*(const float mul)
{
	for (int i = 0; i < this->m_rows; i++)
	{
		for (int j = 0; j < this->m_column; j++)
		{
			this->m_data[i][j] *= mul;
		}
	}
	return *this;
}

const Matrix& Matrix::operator*(const float mul) const
{
	for (int i = 0; i < this->m_rows; i++)
	{
		for (int j = 0; j < this->m_column; j++)
		{
			this->m_data[i][j] *= mul;
		}
	}
	return *this;
}

Matrix& Matrix::operator*=(const float mul)
{
	for (int i = 0; i < this->m_rows; i++)
	{
		for (int j = 0; j < this->m_column; j++)
		{
			this->m_data[i][j] *= mul;
		}
	}
	return *this;
}

const Matrix& Matrix::operator*=(const float mul) const
{
	for (int i = 0; i < this->m_rows; i++)
	{
		for (int j = 0; j < this->m_column; j++)
		{
			this->m_data[i][j] *= mul;
		}
	}
	return *this;
}

void Matrix::SetSubMatrix(const Matrix &mat, int fromRow, int fromCol)
{
	assert(this->m_column >= mat.ColumnSize() && this->m_rows >= mat.RowSize());
	for (int i = fromRow,k=0; i < this->m_rows && k<mat.RowSize(); i++,k++)
	{
		for (int j = fromCol,l=0; j < this->m_column && l<mat.ColumnSize(); j++,l++)
		{
			this->m_data[i][j] = mat[k][l];
		}
	}
}

void Matrix::SetColumnOfMatrixFromRow(const Matrix &mat, int col, int row)
{
	for (int i = row, j = 0; i < mat.RowSize(); i++, j++)
	{
		float a = mat[i][col];;
		this->m_data[j][0] = mat[i][col];
	}
}

Matrix Matrix::operator*(const Matrix &rhs)
{
	Matrix ret = this->Multiply(rhs);
	return ret;
}
