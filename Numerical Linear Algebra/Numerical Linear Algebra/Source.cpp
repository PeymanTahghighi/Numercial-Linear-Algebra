#include<iostream>
#include"GaussElimination.h"

int main()
{
	
	Matrix b(3,1);
	b.UpdateData(0, 0, 8);
	b.UpdateData(1, 0, -1);
	b.UpdateData(2, 0, 13);

	std::cout << "B:\n" << b.MatrixDataToString();
	int a = -1;
	std::cout << std::abs(a);
	Matrix m(3, 3);
	m.UpdateData(0, 0, 1);
	m.UpdateData(0, 1, -1);
	m.UpdateData(0, 2, 3);
	m.UpdateData(1, 0, 2);
	m.UpdateData(1, 1, 1);
	m.UpdateData(1, 2, -1);
	m.UpdateData(2, 0, 3);
	m.UpdateData(2, 1,-6);
	m.UpdateData(2, 2,2);
	std::cout << "\n\nM:\n" << m.Inverse().MatrixDataToString();
	//GaussElimination ge(m);
	//ge.LUFactorization();

	//std::cout<<m.Inverse().MatrixDataToString();
	GaussElimination ge(m, b);
	//ge.SetCompletePivoting(true);
	//ge.SetScalling(true);
	//ge.Calculate();
	
	getchar();
	return 0;
}