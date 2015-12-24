#include "IterativeMethods.h"
#include<iostream>
#include<vector>


IterativeMethods::IterativeMethods(Matrix &m,Matrix &b)
{
	this->m_matrix = m;
	this->m_b = b;
}


IterativeMethods::~IterativeMethods()
{
}

void IterativeMethods::CalculateUsingJacobiAlgorithm()
{
	float * answers = new float[this->m_b.RowSize()];
	std::vector<float*> iterativeLevels;
	iterativeLevels.resize(this->m_iteration+1);
	iterativeLevels[0] = answers;
	iterativeLevels.at(0)[0] = 1;
	for (int i = 0; i < this->m_b.RowSize(); i++)
	{
		answers[i] = 0.0f;
	}
	for (int i = 0; i < this->m_iteration; i++)
	{
		float * nanswers = new float[this->m_b.RowSize()];
		for (int j = 0; j < this->m_b.RowSize(); j++)
		{
			float sum = 0.0f;
			for (int k = 0; k < this->m_b.RowSize(); k++)
			{
				if (k == j) continue;
				//std::cout << "\nIterative Level:" << i << " Answer:  " << k << " Is  " << iterativeLevels.at(i)[k];
				sum += this->m_matrix[j][k] * iterativeLevels.at(i)[k];
			}
			nanswers[j] = (1 / this->m_matrix[j][j])*((this->m_b[j][0]) - sum);
		}

	
		iterativeLevels.at(i+1) = nanswers;
		std::cout << "\n\n*======= Answers At Level " << i << " Is:";
		for (int j = 0; j < this->m_b.RowSize(); j++)
		{
			std::cout << "\n Answer " << j << " : " << nanswers[j];
		}
		std::cout << "\n=======*";
		
	}
	std::cout << "\n\n*================Final Answers:\n";
	for (int i = 0; i < this->m_matrix.RowSize(); i++)
	{
		std::cout << "x" << std::to_string(i) << "=" << iterativeLevels[iterativeLevels.size()-1][i] << "\n";
	}

}

void IterativeMethods::CalculateUsingGaussSidelAlgorithm()
{
	float * answers = new float[this->m_b.RowSize()];
	std::vector<float*> iterativeLevels;
	iterativeLevels.resize(this->m_iteration + 1);
	iterativeLevels[0] = answers;
	iterativeLevels.at(0)[0] = 1;
	for (int i = 0; i < this->m_b.RowSize(); i++)
	{
		answers[i] = 0.0f;
	}
	for (int j = 0; j < iterativeLevels.size();j++)
	{
		iterativeLevels[j] = new float[this->m_b.RowSize()];
		for (int i = 0; i < this->m_b.RowSize(); i++)
		{
			iterativeLevels[j][i] = 0.0f;
		}
	}
	for (int i = 1; i < this->m_iteration+1; i++)
	{
		float * nanswers = new float[this->m_b.RowSize()];

		for (int l = 0; l< this->m_b.RowSize(); l++) nanswers[l] = 0.0f;

		for (int j = 0; j < this->m_b.RowSize(); j++)
		{
			float sum = 0.0f;
			for (int k = 0; k < this->m_b.RowSize(); k++)
			{
				if (k == j) continue;

				if (k < j)
				{
					float a = nanswers[k];
					sum += this->m_matrix[j][k] * nanswers[k];
				}
				else
				{
					float a = iterativeLevels.at(i - 1)[k];
					sum += this->m_matrix[j][k] * iterativeLevels.at(i - 1)[k];
				}
				/*for (int l = j; l >= 0; l--)
				{

				}

				if (j == 0)
				{
					if (i != 0)
					{
						float a = iterativeLevels.at(i - 1)[k];
						sum += this->m_matrix[j][k] * iterativeLevels.at(i - 1)[k];
					}
						
					else
						sum += this->m_matrix[j][k] * 0.0f;
				}
					
				else
					sum += this->m_matrix[j][k] * nanswers[k];*/
			}
			float a = this->m_b[j][0];
			float b = this->m_matrix[j][j];
			nanswers[j] = (1 / this->m_matrix[j][j])*((this->m_b[j][0]) - sum);
		}


		iterativeLevels.at(i) = nanswers;
		std::cout << "\n\n*======= Answers At Level " << i << " Is:";
		for (int j = 0; j < this->m_b.RowSize(); j++)
		{
			std::cout << "\n Answer " << j << " : " << nanswers[j];
		}
		std::cout << "\n=======*";

	}
	std::cout << "\n\n*================Final Answers:\n";
	for (int i = 0; i < this->m_matrix.RowSize(); i++)
	{
		std::cout << "x" << std::to_string(i) << "=" << iterativeLevels[iterativeLevels.size() - 1][i] << "\n";
	}
}
