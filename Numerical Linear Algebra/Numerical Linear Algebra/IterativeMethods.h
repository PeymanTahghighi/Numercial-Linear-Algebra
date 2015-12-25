#pragma once
#include"Matrix.h"

class IterativeMethods
{
public:

	IterativeMethods(Matrix &m,Matrix &b);
	IterativeMethods() = delete;
	~IterativeMethods();

public:

	void SetIteration(int it)			{ this->m_iteration = it; }
	int GetIteration()					{ return this->m_iteration; }

	void SetTelorance(float t)			{ this->m_telorance = t; }
	float GetTelorance()				{ return this->m_telorance; }

	void SetLog(bool b)					{ this->m_log = b; }
	bool GetLog()						{ return this->m_log; }

public:
	std::vector<float> CalculateUsingJacobiAlgorithm();
	std::vector<float> CalculateUsingGaussSidelAlgorithm();



private:
	Matrix m_matrix, m_b;
	int m_iteration;
	float m_telorance;
	bool m_log;
};

