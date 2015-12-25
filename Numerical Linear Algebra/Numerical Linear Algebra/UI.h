#pragma once
#include"GaussElimination.h"
class UI
{
public:
	UI();
	~UI();

public:
	void Initialize();
	void ClearScreen();

	void RenderMainMenu();
	void ProcessInput(int i);

private:
	Matrix m_matrix;
	Matrix m_b;
};

