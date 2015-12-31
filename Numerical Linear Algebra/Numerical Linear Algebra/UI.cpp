#include "UI.h"
#include<iostream>
#include"IterativeMethods.h"


UI::UI()
{
}


UI::~UI()
{
}

void UI::Initialize()
{
	ClearScreen();
	std::cout << "------------------------------------------------------";
	std::cout << "\nEnter Number of Rows:";
	int r;
	std::cin >> r;
	std::cout << "\nEnter Number of Columns:";
	int c;
	std::cin >> c;
	
	Matrix m(r, c);
	this->m_matrix = m;
	ClearScreen();
	float v = 0.0f;
	std::cout << "------------------------------------------------------";
	std::cout << "\nEnter Matrix M Values\n\n";
	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < c; j++)
		{
			std::cout << "Enter Element M[" << i << "][" << j << "]:";
			std::cin >> v;
			this->m_matrix[i][j] = v;
		}
	}

	ClearScreen();
	Matrix b(r, 1);
	this->m_b = b;
	std::cout << "------------------------------------------------------";
	std::cout << "\nEnter Matrix B Values\n\n";

	for (int i = 0; i < r; i++)
	{
		std::cout << "Enter Element B[" << i << "][" << 0 << "]:";
		std::cin >> v;
		this->m_b[i][0] = v;
	}
	RenderMainMenu();
}

void UI::ClearScreen()
{
	system("CLS");
}

void UI::RenderMainMenu()
{
	ClearScreen();
	std::cout << "------------------------------------------------------";
	std::cout << "\nMatrix M:\n" << this->m_matrix.MatrixDataToString();
	std::cout << "\nMatrix B:\n" << this->m_b.MatrixDataToString();
	std::cout << "\n\n*========== Choices ====================*\n";
	std::cout << "1-Calculate System Using Gauss Elemination.";
	std::cout << "\n2-Calculate Invert Of Matrix M.";
	std::cout << "\n3-Calculate Upper Triangular Of Matrix M.";
	std::cout << "\n4-Calculate Determinante Of Matrix M.";
	std::cout << "\n5-Calculate QR Factorization Of Matrix M.";
	std::cout << "\n6-Calculate Characteristic Polynomial Of M.";
	std::cout << "\n7-Calculate System Using Jacobi Method.";
	std::cout << "\n8-Calculate System Using Gauss-Sidel Method.";
	std::cout << "\n9-Calculate Least Square Answers.";
	std::cout << "\n10-Enter Eigenvalues To Create Matrix.(PROJECT)";
	std::cout << "\n0-To Exit.";
	int i;
	std::cout << "\n\nChoice:";
	std::cin >> i;
	ProcessInput(i);
}

void UI::ProcessInput(int i)
{
	switch (i)
	{
	case 1:
	{
		ClearScreen();
		bool b;
		GaussElimination ge(this->m_matrix, this->m_b);
		std::cout << "------------------------------------------";
		std::cout << "\nSet Scalling:";
		std::cin >> b;
		ge.SetScalling(b);
		std::cout << "\nSet Complete Pivoting:";
		std::cin >> b;
		ge.SetCompletePivoting(b);
		if (!b)
		{
			std::cout << "Set Partial Pivoting:";
			std::cin >> b;
			ge.SetPartialPivoting(b);
		}
		ClearScreen();
		
		std::vector<float> ans = ge.Calculate();
		for (int i = 0; i < ans.size(); i++)
		{
			std::cout << "x" << std::to_string(i) << "=" << ans[i] << "\n";
		}
		getchar();
		getchar();
		RenderMainMenu();
	}
		break;
	case 2:
	{
		ClearScreen();
		if (this->m_matrix.IsSingular())
		{

			std::cout << "ERROR : Matrix Is Singular!";
			getchar();
			getchar();
			RenderMainMenu();
			break;
		}
		std::cout << "\n\n==*Matrix Inverse\n";
		std::cout << this->m_matrix.Inverse().MatrixDataToString();
		getchar();
		getchar();
		RenderMainMenu();
	}
	break;
	case 3:
	{
		GaussElimination ge(this->m_matrix);
		ge.CalulateUpper();
		ClearScreen();
		std::cout << "--------------------------------\n";
		std::cout << "Upper Triangular Matrix M Is:\n\n";
		std::cout << ge.GetUpperMatrix().MatrixDataToString();
		getchar();
		getchar();
		RenderMainMenu();
	}
	break;
	case 4:
	{
		ClearScreen();
		std::cout << "---------------------------------------";
		std::cout << "\nDeteminante Of Matrix M Is:" << this->m_matrix.Determinate();
		getchar();
		getchar();
		RenderMainMenu();
	}
	break;
	case 5:
	{
		ClearScreen(); 
		std::cout << "---------------------------------------\n";
		this->m_matrix.QRFactorization(nullptr,nullptr);
		getchar();
		getchar();
		RenderMainMenu();
	}
	break;
	case 6:
	{
		ClearScreen();
		std::cout << "---------------------------------------\n";
		this->m_matrix.GetCharacteristicPolynomial();
		getchar();
		getchar();
		RenderMainMenu();
	}
	break;

	case 7:
	{
		ClearScreen();
		int i = 0;
		std::cout << "Number Of Iterations:";
		std::cin >> i;
		std::cout << "---------------------------------------\n";
		IterativeMethods it(this->m_matrix, this->m_b);
		it.SetIteration(i);
		it.SetLog(true);
		it.CalculateUsingJacobiAlgorithm();
		getchar();
		getchar();
		RenderMainMenu();
	}
	break;

	case 8:
	{
		ClearScreen();
		int i = 0;
		std::cout << "Number Of Iterations:";
		std::cin >> i;
		std::cout << "---------------------------------------\n";
		IterativeMethods it(this->m_matrix, this->m_b);
		it.SetIteration(i);
		it.SetLog(true);
		it.CalculateUsingGaussSidelAlgorithm();
		getchar();
		getchar();
		RenderMainMenu();
	}
	break;

	case 9:
	{
		ClearScreen();
		std::cout << "------------------------------------------";
		if (this->m_matrix.RowSize() == this->m_matrix.ColumnSize())
		{
			std::cout << "ERROR : Cannot Find Least Squares Answers!";
			getchar();
			getchar();
			RenderMainMenu();
		}
		GaussElimination ge(this->m_matrix, this->m_b);
		ge.CalculateLeastSquares();
		getchar();
		getchar();
		RenderMainMenu();
	}
	break;

	case 10:
	{
		ClearScreen();
		std::cout << "------------------------------";
		std::cout << "\nEnter Number Of EigenValues";
		int i;
		std::cin >> i;
		std::vector<float> eigen;
		int n = 0;
		while (n < i)
		{
			float number;
			std::cout << "Enter Number " << n << " :";
			std::cin >> number;
			eigen.push_back(number);
			n++;
		}
		Matrix m(i, i);
		m.SetDiagonalElements(eigen);
		ClearScreen();
		
		std::vector<float> val;
		val = m.GetCharacteristicPolynomial();
		std::reverse(val.begin(), val.end());
		for (auto & member : val)
			member *= -1;
		std::cout << "\n\nNow We Calculate The Matrix\n\nPress Enter To Continue...";
		getchar();
		getchar();
		Matrix m1(i, i);
		m1.SetUnderRowElements(1.0f);
		m1.SetRow(0, val);
		ClearScreen();
		std::cout << "---------------------------------------------";
		std::cout << "\nMatrix That Have Enter Eigenvalues is:\n\n";
		std::cout << m1.MatrixDataToString();
		std::cout << "\n\nChecking Characteristic Polynomial Of New Matrix:\n";
		m1.GetCharacteristicPolynomial();
		std::cout << "Now We Use Deflation To Calculate All EigenValues\n\nPress Enter To Continue...";
		getchar();
		ClearScreen();
		std::cout << "---------------------------------------";
		std::cout << "\n\nEnter Precision:";
		float p;
		std::cin >> p;
		Matrix e(i, 1);
		std::vector<float> ans =  m1.GetEigenValuesUsingDeflation(p);
		ClearScreen();
		std::cout << "**=========================**";
		for (auto const & member : ans)
		{
			std::cout << "\n\n==*\n EigenValue : " << member;
			std::cout << "\n*==";
		}
		getchar();
		getchar();
		RenderMainMenu();
	}
	break;
	case 1100:
	{
		ClearScreen();
		GaussElimination ge(this->m_matrix);
		ge.LUFactorization();
		getchar();
		getchar();
		RenderMainMenu();
	}
	break;
	case 0:
		break;
	}
}