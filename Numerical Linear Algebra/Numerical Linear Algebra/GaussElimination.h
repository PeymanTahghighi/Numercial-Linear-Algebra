#pragma once
#include"Matrix.h"
#include<vector>

class GaussElimination
{
public:
	struct Operation
	{
		float m;
		int row;
		int addedToRow;
	};
	struct SwapRecord
	{
		int swap;
		int swapBy;
	};

	
public:
	GaussElimination(Matrix m,Matrix b);
	GaussElimination(Matrix m);
	~GaussElimination();

public:
	void SetScalling(bool b) { this->m_scalling = b; }
	bool GetScalling() const { return this->m_scalling; }

	void SetPartialPivoting(bool b) { this->m_partialPivoting = b; }
	bool GetPartialPivoting() const { return this->m_partialPivoting; }

	void SetCompletePivoting(bool b) { this->m_completePivoting = b; }
	bool GetCompletePivoting() const { return this->m_completePivoting; }

	Matrix &GetUpperMatrix() { return this->m_matrix; }
	int &GetSignChangeCounter()  { return this->m_signChange; }

public:
	
	std::vector<float> Calculate();
	void CalulateUpper();
	void LUFactorization();
	void CalculateLeastSquares();

private:
	void DoScalling();

private:
	Matrix m_matrix, m_b;
	std::vector<Operation> m_operations;
	int m_signChange;

	bool m_scalling, m_partialPivoting, m_completePivoting;
};

