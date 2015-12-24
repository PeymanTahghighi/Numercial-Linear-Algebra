#include<iostream>
#include"GaussElimination.h"
#include"IterativeMethods.h"

int main()
{
	Matrix B(30, 1);
	B[0][0] = 8;
	B[1][0] = -14;
	B[2][0] = 27;


	Matrix A(2, 2);
	A[0][0] = 2;
	A[0][1] = -12;
	A[1][0] = 1;
	A[1][1] = -5;



	/*A[1][0] =56;
	A[1][1] =45;
	A[1][2] = 89;
	A[1][3] = 6;


	A[2][0] = 87;
	A[2][1] =65;
	A[2][2] = 12;
	A[2][3] = 6;

	A[3][0] = 6;
	A[3][1] = 6;
	A[3][2] = 6;
	A[3][3] = 6;*/



	A.GetCharacteristicPolynomial();

	float max = A.GetEigenvalues(0.00001);
	std::cout << "/n/n**==" << max;
	
	getchar();
	return 0;
}