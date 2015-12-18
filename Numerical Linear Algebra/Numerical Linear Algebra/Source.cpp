#include<iostream>
#include"GaussElimination.h"

int main()
{
	Matrix B(4, 1);
	B[0][0] = 2;
	B[1][0] = 4;
	B[2][0] = 0;
	B[3][0] = 1;

	Matrix A(3, 2);
	A[0][0] = 1;
	A[0][1] = 1;

	A[1][0] =1;
	A[1][1] = 2;

	A[2][0] = 1;
	A[2][1] = 3;

	/*A[2][0] = 2;
	A[2][1] = 0;
	A[2][2] = 1;*/

	/*A[3][0] = 0.9134;
	A[3][1] = 0.9575;
	A[3][2] = 0.4854;

	A[4][0] = 0.6324;
	A[4][1] = 0.9649;
	A[4][2] = 0.8003;*/


	std::cout << "\n\n\n" << A.MatrixDataToString();;
	//std::cout << "\n\n\n" << A.GetColumnInMatrix(1).MatrixDataToString();;

A.QRFactorization();;
	//GaussElimination ge(m);
	//ge.LUFactorization();

	//std::cout<<m.Inverse().MatrixDataToString();
	/*GaussElimination ge(q, b);
	ge.SetCompletePivoting(true);
	ge.SetScalling(true);
	ge.Calculate();*/
	
	getchar();
	return 0;
}