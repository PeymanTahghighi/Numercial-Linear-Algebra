#include "GaussElimination.h"
#include<iostream>
#include<string>

GaussElimination::GaussElimination(Matrix m,Matrix b)
{
	this->m_matrix = m;
	this->m_b = b;
}

GaussElimination::GaussElimination(Matrix m)
{
	this->m_matrix = m;
}

GaussElimination::~GaussElimination()
{
}

std::vector<float> GaussElimination::Calculate()
{ 
	assert(this->m_matrix.RowSize() == this->m_b.RowSize());
	if (this->m_scalling)
		DoScalling();
	int k = 1;
	int maxi = 0;
	int maxj = 0;
	std::vector<SwapRecord> swaps;
	for (int i = 0; i < this->m_matrix.RowSize(); i++)
	{
		if (m_partialPivoting)
		{
			int m = this->m_matrix.FindMaxIndiceInColumn(i,i);
			this->m_matrix.SwapRows(m, i);
			this->m_b.SwapRows(m,i);
		}
		else if (this->m_completePivoting)
		{
			this->m_matrix.FindMaxInBlock(i, &maxi, &maxj);
			this->m_matrix.SwapRows(maxi, i);
			this->m_b.SwapRows(maxi, i);
			this->m_matrix.SwapColumns(maxj, i);
			SwapRecord sr;
			sr.swap = maxj;
			sr.swapBy = i;
			swaps.push_back(sr);
			
		}

		for (int j = k; j < this->m_matrix.RowSize(); j++)
		{
			Operation o;
			float m = 0.0f;
			m = -(this->m_matrix[j][i]) / (this->m_matrix[i][i]);
			this->m_matrix.MultiplyRowByAndAddToRow(i, m, j);
			this->m_b.MultiplyRowByAndAddToRow(i, m, j);
			
			o.m = m;
			o.row = i;
			o.addedToRow = j;
			this->m_operations.push_back(o);
		}
		k++;
	}
	float *answers = new float[this->m_matrix.RowSize()];
	k = this->m_matrix.RowSize() - 1;
	for (int i = this->m_matrix.RowSize() - 1; i >= 0; i--)
	{
		float sum = 0;
		for (int j = this->m_matrix.ColumnSize()-1; j >= k; j--)
		{
			if (i != j)
			{
				sum += answers[j] * this->m_matrix[i][j];
			}
			else
			{
				answers[j] = (this->m_b[k][this->m_b.ColumnSize() - 1] - sum) / this->m_matrix[i][j];
			}
		}
		k--;
	}

	std::vector<float> ret;
	
	for (int i = swaps.size() - 1; i >= 0;i--)
	{
		float temp = answers[swaps[i].swap];
		answers[swaps[i].swap] = answers[swaps[i].swapBy];
		answers[swaps[i].swapBy] = temp;

	}
	
	for (int i = 0; i < this->m_matrix.RowSize(); i++)
	{
		ret.push_back(answers[i]);
	}

	return ret;
}

void GaussElimination::DoScalling()
{
	for (int i = 0; i < this->m_matrix.RowSize(); i++)
	{
		float max = 0;
		for (int j = 0; j < this->m_matrix.ColumnSize(); j++)
		{
			if (std::abs(this->m_matrix[i][j])>max) max = std::abs(this->m_matrix[i][j]);
		}
		for (int j = 0; j < this->m_b.ColumnSize(); j++)
		{
			if (std::abs(this->m_b[i][j])>max) max = std::abs(this->m_b[i][j]);
		}
		this->m_matrix.DivideEachElementInRowBy(i, max);
		this->m_b.DivideEachElementInRowBy(i, max);
	}
}

void GaussElimination::CalulateUpper()
{
	this->m_matrix.ScaleMatrix();
	int k = 1;
	int maxi = 0;
	int maxj = 0;
	for (int i = 0; i < this->m_matrix.RowSize(); i++)
	{
		
		this->m_matrix.FindMaxInBlock(i, &maxi, &maxj);
		this->m_matrix.SwapRows(maxi, i);
		this->m_matrix.SwapColumns(maxj, i);
		//if we really change two rows.
		//if (m != i)
			this->m_signChange+=2;

		
		for (int j = k; j < this->m_matrix.RowSize(); j++)
		{
			float m = 0.0f;

			m = -(this->m_matrix[j][i]) / (this->m_matrix[i][i]);
			this->m_matrix.MultiplyRowByAndAddToRow(i, m, j);
		}
		k++;
	}
}

void GaussElimination::LUFactorization()
{
	std::vector<Matrix> matricesL;
	std::vector<Matrix> matricesP;
	
	int k = 1;
	int maxi = 0;
	int maxj = 0;
	for (int i = 0; i < this->m_matrix.RowSize(); i++)
	{
		Matrix P(this->m_matrix.RowSize(), this->m_matrix.ColumnSize());
		P.LoadIdentity();
		int m = this->m_matrix.FindMaxIndiceInColumn(i, i);
		this->m_matrix.SwapRows(m, i);
		P.SwapRows(m, i);
		matricesP.push_back(P);

		Matrix M(this->m_matrix.RowSize(), this->m_matrix.ColumnSize());
		M.LoadIdentity();
		for (int j = k; j < this->m_matrix.RowSize(); j++)
		{
			float m = 0.0f;
			
			m = -(this->m_matrix[j][i]) / (this->m_matrix[i][i]);
			this->m_matrix.MultiplyRowByAndAddToRow(i, m, j);
			M.MultiplyRowByAndAddToRow(i, (m)*-1, j);
		}

		matricesL.push_back(M);
		std::cout << "\n\n Permutation\n" << P.MatrixDataToString();
		std::cout << "\n\nM" << M.MatrixDataToString() << "\n\n";
		std::cout << "\n\nMatrix" << this->m_matrix.MatrixDataToString() << "\n\n";
		k++;
	}

	std::cout << "\n\n*===Matrix U\n";
	std::cout << this->m_matrix.MatrixDataToString();
	std::cout << "\n===*\n\n";

	Matrix Multiply(this->m_matrix.RowSize(), this->m_matrix.ColumnSize());
	Multiply = matricesL[0];
	for (unsigned int i = 1; i <matricesL.size(); i++)
	{
		Multiply = Multiply.Multiply(matricesL[i]);
	}

	Matrix permu(this->m_matrix.RowSize(), this->m_matrix.ColumnSize());
	permu = matricesP[0];
	for (unsigned int i = 1; i <matricesP.size(); i++)
	{
		permu = permu.Multiply(matricesP[i]);
	}

	std::cout << "\n\n*===MATRIX L\n\n";
	std::cout <<Multiply.MatrixDataToString();
	std::cout << "\n\n===*\n\n";

	std::cout << "\n\n*===MATRIX P\n\n";
	std::cout << permu.MatrixDataToString();
	std::cout << "\n\n===*\n\n";

	std::cout << "\n\n==*Cheking Answer:\n\n";
	Matrix ans = Multiply.Multiply(this->m_matrix);
	ans = permu*ans;
	std::cout << ans.MatrixDataToString();
	std::cout << "\n\n===*\n\n";

}

void GaussElimination::CalculateLeastSquares()
{
	Matrix aTa(this->m_matrix.ColumnSize(), this->m_matrix.ColumnSize());
	Matrix aTb(this->m_matrix.RowSize(), this->m_b.ColumnSize());
	Matrix T;
	T = this->m_matrix.Transpose();;
	//std::cout << "\n\n*** matrix is:\n" << this->m_matrix.MatrixDataToString();
	//std::cout << "\n\n*** matrix is:\n" << T.MatrixDataToString();
	aTa = T.Multiply(this->m_matrix);
	//std::cout << "\n\n*** ATA is:\n" << aTa.MatrixDataToString();
	aTb = T.Multiply(this->m_b);
	//std::cout << "\n\n*** ATB is:\n" << aTb.MatrixDataToString();
	GaussElimination ge(aTa, aTb);
	ge.DoScalling();
	ge.SetCompletePivoting(true);
	std::vector<float> ans= ge.Calculate();
	std::cout << "\nLeast Square Answers Are:\n";
	for (int i = 0; i < ans.size(); i++)
	{
		std::cout << "x" << std::to_string(i) << "=" << ans[i] << "\n";
	}
}