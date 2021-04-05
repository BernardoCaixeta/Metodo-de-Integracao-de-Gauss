// ConsoleApplication5.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Gauss.h"
#include <string>
#include <direct.h>
#include "iostream" 
#include "fstream"
using std::ofstream;
using std::ifstream;
using namespace std;


Gauss::Gauss()
{

	grau = 0;
	NumPto = 0;
	PtoInteg = 0;
	PesoGauss = 0;
	X1 = 0;
	X2 = 0;
	X3 = 0;
	Y1 = 0;
	Y2=  0;
	Y3 = 0;
	Matriz_Jacob[2][2] = 0;
	Matriz_Jacob_Inv[2][2] = 0;
	Determinante = 0;
	Resultado = 0;
	

}

// Pergunta o grau do polinômio ao usuário.
void Gauss::NumMaxGrauParPtoGauss()
{
	double n;
	
	cout << "Programa desenvolvido para avaliacao da diciplina 'Introducao ao c++' minitrada pelo professor Eduardo Carmo no PLE2020" << "\n";
	
	while(true) {

		cout << "Digite o grau do polinomio " << "\n";
		
		cin >> n;
		

		if (0<n && n<7)
		{
			cout << "O grau para integracao de Gauss e:";
			cout << n << "\n";
			Setgrau(n); // Setter do grau do polinômio.
			break;
		}
		else
		{
			cout << "Valor invalido para o grau do polinomio" << "\n";
		}

	}
}

void Gauss::TesteDiretorio() {
	
	string    DiscoValido, DiretorioValido, Nome, Arquivo, Endereco;
	ofstream  arquivo;
	double X, Y, Peso;
	double L0, L1, L2;
	double Vert_X[3], Vert_Y[3];

	X = 0;
	Y = 0;
	Peso = 0;
	Vert_X[0] = GetX1();
	Vert_X[1] = GetX2();
	Vert_X[2] = GetX3();
	Vert_Y[0] = GetY1();
	Vert_Y[1] = GetY2();
	Vert_Y[2] = GetY3();
	
	
	while(true) {
		cout << "Digite o disco no qual deseja salvar os dados do programa:" << "\n";
		cin >> DiscoValido;
		cout << "Digite a pasta na qual deseja salvar os dados do programa:" << "\n";
		cin >> DiretorioValido;
		cout << "Digite o nome do arquivo no qual deseja salvar os dados do programa:" << "\n";
		cin >> Nome;
		
		Endereco = DiscoValido + ":/" + DiretorioValido;
		Arquivo = DiscoValido + ":/" + DiretorioValido + "/" + Nome + ".txt";
		
		_mkdir(Endereco.c_str());
		arquivo.open(Arquivo.c_str());
		
		if (arquivo.fail())
		{
			cout << "O disco e invalido" << "\n";
			cout << "Tente novamente." << "\n";
		}
		else{
			break;
		}
	}

	cout << "\n" << "Os dados do programa serão salvos em: " <<  Arquivo << "\n";

	arquivo << "Pontos: " << "\n";
	arquivo << X1 << " " << Y1 << "\n";
	arquivo << X2 << " " << Y2 << "\n";
	arquivo << X3 << " " << Y3 << "\n" << "\n";

	arquivo << "Matriz Jacobiana: " << "\n";
	arquivo << Matriz_Jacob[0][0] << " " << Matriz_Jacob[0][1] << "\n";
	arquivo << Matriz_Jacob[1][0] << " " << Matriz_Jacob[1][1] << "\n" << "\n";

	arquivo << "Matriz Jacobiana Invesa: " << "\n";
	arquivo << Matriz_Jacob_Inv[0][0] << " " << Matriz_Jacob_Inv[0][1] << "\n";
	arquivo << Matriz_Jacob_Inv[1][0] << " " << Matriz_Jacob_Inv[1][1] << "\n" << "\n";

	arquivo << "Determinante da matriz: " << Determinante << "\n" << "\n";

	arquivo << "Número total de pontos de Gauss: " << NumPto << "\n" << "\n";

	for (int i = 0; i < NumPto; i++) {

		Peso = GetPesoGauss()[i];
		L0 = GetPtoInteg()[i][0];
		L1 = GetPtoInteg()[i][1];
		L2 = GetPtoInteg()[i][2];
		X = Leva_L0_L1_L2_a_X(L0, L1, L2, Vert_X);
		Y = Leva_L0_L1_L2_a_Y(L0, L1, L2, Vert_Y);
		arquivo << " Pontos em coordenadas globais:     " << X << "  " << Y << "\n";
		arquivo << " Pontos em coordenadas área:    " << L0 << "  " << L1 << "  " << L2 << "\n";
		arquivo << " Peso para integração:     " << Peso <<  "\n";
		arquivo << "\n";
		
	}
	
	arquivo << "\n";
	arquivo << "Resultado da integral:  " << "  " << Resultado;

	arquivo.close();

	DalocPtoInteg(Var);
	DalocPesoGauss();

}

void Gauss::NumPtoInteg() 
{	
	int  N1, N2;
	double Pto, Peso;
	
	
	double  CordLocal_PtoDeGauss_Grau_1[1][3] = { 0.333333333333333, 0.333333333333333, 0.333333333333333 };
	double  CordLocal_PtoDeGauss_Grau_2[3][3] = { { 0.5, 0.5, 0.5 }, { 0, 0.5, 0.5 }, { 0.5, 0, 0.5 } };
	double  CordLocal_PtoDeGauss_Grau_3[4][3] = { { 0.333333333333333, 0.333333333333333, 0.333333333333333 }, { 0.6, 0.2, 0.2 }, { 0.2, 0.6, 0.2 }, { 0.2, 0.2, 0.6 } };
	double  CordLocal_PtoDeGauss_Grau_4[6][3] = { { 0.445948490915965, 0.445948490915965, 0.10810301816807 }, { 0.445948490915965, 0.10810301816807, 0.445948490915965 },
	{ 0.10810301816807, 0.445948490915965, 0.445948490915965 }, { 0.091576213509771, 0.091576213509771, 0.816847572980458 },
	{ 0.091576213509771, 0.816847572980458, 0.091576213509771 }, { 0.816847572980458, 0.091576213509771, 0.091576213509771 } };
	double  CordLocal_PtoDeGauss_Grau_5[7][3] = { { 0.333333333333333, 0.333333333333333, 0.333333333333333 }, { 0.101286507323456, 0.101286507323456, 0.797426985353087 },
	{ 0.101286507323456, 0.797426985353087, 0.101286507323456 }, { 0.797426985353087, 0.101286507323456, 0.101286507323456 },
	{ 0.470142064105115, 0.470142064105115, 0.0597158717897699 }, { 0.470142064105115, 0.0597158717897699, 0.470142064105115 },
	{ 0.0597158717897699, 0.470142064105115, 0.470142064105115 } };
	double  CordLocal_PtoDeGauss_Grau_6[12][3] = { { 0.310352451033785, 0.053145049844816, 0.636502499121399 }, { 0.310352451033785, 0.636502499121399, 0.053145049844816 },
	{ 0.636502499121399, 0.310352451033785, 0.053145049844816 }, { 0.053145049844816, 0.310352451033785, 0.636502499121399 },
	{ 0.053145049844816, 0.636502499121399, 0.310352451033785 }, { 0.636502499121399, 0.053145049844816, 0.310352451033785 },
	{ 0.063089014491502, 0.063089014491502, 0.873821971016996 }, { 0.063089014491502, 0.873821971016996, 0.063089014491502 },
	{ 0.873821971016996, 0.063089014491502, 0.063089014491502 }, { 0.24928674517091, 0.24928674517091, 0.50142650965818 },
	{ 0.24928674517091, 0.50142650965818, 0.24928674517091 }, { 0.50142650965818, 0.24928674517091, 0.24928674517091 } };

	double Peso_PtoDeGauss_Grau_1[1] = {1};
	double Peso_PtoDeGauss_Grau_2[3] = { 0.333333333333333, 0.333333333333333, 0.333333333333333 };
	double Peso_PtoDeGauss_Grau_3[4] = { -0.5625, 0.520833333333333, 0.520833333333333, 0.520833333333333 };
	double Peso_PtoDeGauss_Grau_4[6] = { 0.22338158967801, 0.22338158967801, 0.22338158967801, 0.109951743655322, 0.109951743655322, 0.109951743655322 };
	double Peso_PtoDeGauss_Grau_5[7] = { 0.225, 0.125939180544827, 0.125939180544827, 0.125939180544827, 0.132394152788506, 0.132394152788506, 0.132394152788506 };
	double Peso_PtoDeGauss_Grau_6[12] = { 0.082851075618374, 0.082851075618374, 0.082851075618374, 0.082851075618374, 0.082851075618374 ,0.082851075618374, 0.050844906370206,
									  0.050844906370206, 0.050844906370206, 0.116786275726378, 0.116786275726378, 0.116786275726378 };

	

	switch (grau)
	{ 
	
	case 1:    { 
		N1 = 1;
		N2 = 3;
		SetVar(N1);
		SetNumPto(N1);
		
		AlocPtoInteg(N1, N2);
		AlocPesoGauss(N1);

		for (int i = 0; i < N1; i++)
		{
			for (int j = 0; j < N2; j++)
			{
				Pto = CordLocal_PtoDeGauss_Grau_1[i][j];
				SetPtoInteg(Pto, i, j);		
			}
		};

		for (int i = 0; i < N1; i++)
		{
			Peso = Peso_PtoDeGauss_Grau_1[i];
			SetPesoGauss(Peso, i);
		};
		
			break;    } //Fim do escopo
	case 2:  { //Início do escopo
		N1 = 3;
		N2 = 3;
		SetVar(N1);
		SetNumPto(N1);

		AlocPtoInteg(N1, N2);
		AlocPesoGauss(N1);

		for (int i = 0; i < N1; i++)
		{
			for (int j = 0; j < N2; j++)
			{
				Pto = CordLocal_PtoDeGauss_Grau_2[i][j];
				SetPtoInteg(Pto, i, j);
			}
		};

		for (int i = 0; i < N1; i++)
		{
			Peso = Peso_PtoDeGauss_Grau_2[i];
			SetPesoGauss(Peso, i);
		};
		
			break;  } //Fim do escopo

	case 3:  { //Início do escopo
		N1 = 4;
		N2 = 3;
		SetVar(N1);
		SetNumPto(N1);

		AlocPtoInteg(N1, N2);
		AlocPesoGauss(N1);

		for (int i = 0; i < N1; i++)
		{
			for (int j = 0; j < N2; j++)
			{
				Pto = CordLocal_PtoDeGauss_Grau_3[i][j];
				SetPtoInteg(Pto, i, j);
			}
		};

		for (int i = 0; i < N1; i++)
		{
			Peso = Peso_PtoDeGauss_Grau_3[i];
			SetPesoGauss(Peso, i);
		};

		break;   }

	case 4:  { //Início do escopo
		N1 = 6;
		N2 = 3;
		SetVar(N1);
		SetNumPto(N1);

		AlocPtoInteg(N1, N2);
		AlocPesoGauss(N1);

		for (int i = 0; i < N1; i++)
		{
			for (int j = 0; j < N2; j++)
			{
				Pto = CordLocal_PtoDeGauss_Grau_4[i][j];
				SetPtoInteg(Pto, i, j);
			}
		};

		for (int i = 0; i < N1; i++)
		{
			Peso = Peso_PtoDeGauss_Grau_4[i];
			SetPesoGauss(Peso, i);
		};

		break;   }

	case 5:  { //Início do escopo
		N1 = 7;
		N2 = 3;
		SetVar(N1);
		SetNumPto(N1);

		AlocPtoInteg(N1, N2);
		AlocPesoGauss(N1);

		for (int i = 0; i < N1; i++)
		{
			for (int j = 0; j < N2; j++)
			{
				Pto = CordLocal_PtoDeGauss_Grau_5[i][j];
				SetPtoInteg(Pto, i, j);
			}
		};

		for (int i = 0; i < N1; i++)
		{
			Peso = Peso_PtoDeGauss_Grau_5[i];
			SetPesoGauss(Peso, i);
		};

		break;   }

	case 6:  { //Início do escopo
		N1 = 12;
		N2 = 3;
		SetVar(N1);
		SetNumPto(N1);

		AlocPtoInteg(N1, N2);
		AlocPesoGauss(N1);

		for (int i = 0; i < N1; i++)
		{
			for (int j = 0; j < N2; j++)
			{
				Pto = CordLocal_PtoDeGauss_Grau_6[i][j];
				SetPtoInteg(Pto, i, j);
			}
		};

		for (int i = 0; i < N1; i++)
		{
			Peso = Peso_PtoDeGauss_Grau_6[i];
			SetPesoGauss(Peso, i);
		};

		break;   }
	
	default:  { //Início do escopo
		cout << " Valor escolido para o grau e invalido " << "\n";
		break;    } //Fim do escopo
	} //Fim do escope do switch(m)

};

void Gauss::Coordenada_Global_Vertice()
{
	double num;

	num = 0;

	cout << "digite a coordenada global X1 do vertice do triangulo:" << "\n";
	cin >> num;
	SetX1(num);

	cout << "digite a coordenada global Y1 do vertice do triangulo:" << "\n";
	cin >> num;
	SetY1(num);

	cout << "digite a coordenada global X2 do vertice do triangulo:" << "\n";
	cin >> num;
	SetX2(num);

	cout << "digite a coordenada global Y2 do vertice do triangulo:" << "\n";
	cin >> num;
	SetY2(num);

	cout << "digite a coordenada global X3 do vertice do triangulo:" << "\n";
	cin >> num;
	SetX3(num);

	cout << "digite a coordenada global Y3 do vertice do triangulo:" << "\n";
	cin >> num;
	SetY3(num);


};

double Gauss::Leva_L0_L1_L2_a_X(double& L0, double& L1, double& L2, double* Vert_X)
{
	double X;
	X = (Vert_X[0] - Vert_X[2])*L0 + (Vert_X[1] - Vert_X[2])*L1 + Vert_X[2];
	return X;
};

double Gauss::Leva_L0_L1_L2_a_Y(double& L0, double& L1, double& L2, double* Vert_Y)
{
	double Y;
	Y = (Vert_Y[0] - Vert_Y[2])*L0 + (Vert_Y[1] - Vert_Y[2])*L1 + Vert_Y[2];
	return Y;
};

// Calcula a inversa da matriz
void Gauss::Matriz()
{
	double		   X1, X2, X3, Y1, Y2, Y3, k[4], det, j[4];
	

	X1 = GetX1();
	X2 = GetX2();
	X3 = GetX3();
	Y1 = GetY1();
	Y2 = GetY2();
	Y3 = GetY3();

	k[0] = X1 - X3;
	k[1] = Y1 - Y3;
	k[2] = X2 - X3;
	k[3] = Y2 - Y3;

	SetMatriz_Jacob(k[0], 0, 0);
	SetMatriz_Jacob(k[1], 0, 1);
	SetMatriz_Jacob(k[2], 1, 0);
	SetMatriz_Jacob(k[3], 1, 1);

	det = k[0] * k[3] - k[1] * k[2];

	SetDeterminante(det);

	j[0] = k[3] / det;
	j[1] = -k[1] / det;
	j[2] = -k[2] / det;
	j[3] = k[0] / det;

	SetMatriz_Jacob_Inv(j[0], 0, 0);
	SetMatriz_Jacob_Inv(j[1], 0, 1);
	SetMatriz_Jacob_Inv(j[2], 1, 0);
	SetMatriz_Jacob_Inv(j[3], 1, 1);
	
	cout << "pontos" << "\n";
	cout << X1 << " " << Y1 << "\n";
	cout << X2 << " " << Y2 << "\n";
	cout << X3 << " " << Y3 << "\n" << "\n" << "\n";

	cout << "matriz" << "\n";
	cout << k[0] << " " << k[1] << "\n";
	cout << k[2] << " " << k[3] << "\n" << "\n" << "\n";

	cout << "matriz inversa" << "\n";
	cout << j[0] << " " << j[1] << "\n";
	cout << j[2] << " " << j[3] << "\n" << "\n" << "\n";

	cout << "determinante  " << det << "\n" << "\n" << "\n";
	


}

void Gauss::Calcula_Integral_Gauss()
{
	int Num_de_Pontos;
	double Soma_Gauss, X, Y, Func, Peso, Det;
	double L0,L1,L2;
	double Vert_X[3], Vert_Y[3];

	Num_de_Pontos = 0;
	Soma_Gauss = 0;
	X = 0;
	Y = 0;
	Func = 0;
	Det = 0;
	Peso = 0;
	Vert_X[0] = GetX1();
	Vert_X[1] = GetX2();
	Vert_X[2] = GetX3();
	Vert_Y[0] = GetY1();
	Vert_Y[1] = GetY2();
	Vert_Y[2] = GetY3();

	Num_de_Pontos = GetNumPto();
	Det = GetDeterminante();

	for (int i = 0; i < Num_de_Pontos; i++)
	{
		Peso = GetPesoGauss()[i];
		L0 = GetPtoInteg()[i][0];
		L1 = GetPtoInteg()[i][1];
		L2 = GetPtoInteg()[i][2];
		X = Leva_L0_L1_L2_a_X(L0, L1, L2, Vert_X);
		Y = Leva_L0_L1_L2_a_Y(L0, L1, L2, Vert_Y);
		cout << " pontos em coordenadas area    " << L0 << "  " << L1 << "  " << L2 << "\n";
		cout << " pontos em coordenadas globais     " << X << "  " << Y << "\n";
		cout << "\n";
		Func = X*X*X*X*X*X*Y;

		Soma_Gauss = Soma_Gauss + (Func*Peso*Det);
		Soma_Gauss = 0.5*Soma_Gauss; // Adicionado: Divide o resultado por 2

	}

	cout << "Resultado integral:  " << Soma_Gauss << "\n";
	SetResultado(Soma_Gauss);

}

void Gauss::Setgrau(int g)
{
	grau = g;
}

void Gauss::SetNumPto(int num)
{
	NumPto = num;
}

void Gauss::SetVar(int num)
{
	Var = num;
}

void  Gauss::SetPtoInteg(double  Pto, int&     N1, int&   N2)
{
			PtoInteg[N1][N2] =	Pto;
		
}

void Gauss::SetPesoGauss(double Peso, int& N1)
{
	PesoGauss[N1] = Peso;
}

void Gauss::SetX1(double num)
{
	X1 = num;
}

void Gauss::SetX2(double num)
{
	X2 = num;
}

void Gauss::SetX3(double num)
{
	X3 = num;
}

void Gauss::SetY1(double num)
{
	Y1 = num;
}

void Gauss::SetY2(double num)
{
	Y2 = num;
}

void Gauss::SetY3(double num)
{
	Y3 = num;
}

void Gauss::SetMatriz_Jacob(double num, int N1, int N2)
{
	Matriz_Jacob[N1][N2] = num;
}

void Gauss::SetMatriz_Jacob_Inv(double num, int N1, int N2)
{
	Matriz_Jacob_Inv[N1][N2] = num;
}

void Gauss::SetDeterminante(double& det)
{
	Determinante = det;
}

void Gauss::SetResultado(double result)
{
	Resultado = result;
}


int Gauss::Getgrau()
{
	return grau;
}

double* Gauss::GetPesoGauss()
{
	return PesoGauss;
}

double** Gauss::GetPtoInteg()
{
	return PtoInteg;
};

int Gauss::GetNumPto()
{
	return NumPto;
}

int Gauss::GetVar()
{
	return Var;
}

double Gauss::GetX1()
{
	return X1;
}

double Gauss::GetX2()
{
	return X2;
}

double Gauss::GetX3()
{
	return X3;
}

double Gauss::GetY1()
{
	return Y1;
}

double Gauss::GetY2()
{
	return Y2;
}

double Gauss::GetY3()
{
	return Y3;
}

double Gauss::GetMatriz_Jacob(int N1, int N2)
{
	return Matriz_Jacob[N1][N2];
}

double Gauss::GetMatriz_Jacob_Inv(int N1, int N2)
{
	return Matriz_Jacob_Inv[N1][N2];
}

double Gauss::GetDeterminante()
{
	return Determinante;
}


void Gauss::AlocPtoInteg(int&   N1, int&  N2)
{          
	try  {
		PtoInteg = new double*[N1];
	}
	catch (bad_alloc)  {
		
		cout << "Alocacão  Falhou." << "\n";
		exit(EXIT_FAILURE);

	}
	for (int i = 0; i < N1; i++) 
		try   {
		PtoInteg[i] = new double[N2];
		}
	catch (bad_alloc)  {
		cout << "Alocacão  Falhou." << "\n";
		exit(EXIT_FAILURE);
		}
}

void Gauss::AlocPesoGauss(int& N1)
{
		try  {
			PesoGauss = new double[N1];
		}
		catch (bad_alloc)  {

			cout << "Alocacão  Falhou." << "\n";
			exit(EXIT_FAILURE);

		}
};

void Gauss::DalocPtoInteg(int& N1)
{
	{        
		for (int i = 0; i < N1; i++)
		{
			delete[] PtoInteg[i];
		}
		delete[] PtoInteg;
		
	}
}

void Gauss::DalocPesoGauss()
{
	delete[]PesoGauss;
}



