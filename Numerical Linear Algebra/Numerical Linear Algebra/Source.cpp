#include<iostream>
#include"GaussElimination.h"
#include"IterativeMethods.h"
#include"UI.h"

int main()
{
	UI * ui = new UI();
	ui->Initialize();
	getchar();
	delete ui;
	return 0;
}