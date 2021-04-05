#include "stdafx.h"
#include "Gauss.h"
#include <string>
#include <direct.h>
#include "iostream" 
#include "fstream"
using std::ofstream;
using std::ifstream;
using namespace std;

int main()
{
	int t;

	Gauss *trabalho1 = new Gauss();

	trabalho1->NumMaxGrauParPtoGauss();

	trabalho1->NumPtoInteg();

	trabalho1->Coordenada_Global_Vertice();

	trabalho1->Matriz();

	trabalho1->Calcula_Integral_Gauss();

	trabalho1->TesteDiretorio();
	
	cout << "digite um numero para encerrar o programa";
	cin >> t;


	return 0;
}
