#ifndef  _GAUSS_H_included
#define  _GAUSS_H_included


class Gauss{
public:
	Gauss();
	void NumMaxGrauParPtoGauss();
	void TesteDiretorio();
	void NumPtoInteg();
	void Coordenada_Global_Vertice();
	void Matriz();
	void Calcula_Integral_Gauss();
	double Leva_L0_L1_L2_a_X(double& L0, double& L1, double& L2, double* Vert_X);
	double Leva_L0_L1_L2_a_Y(double& L0, double& L1, double& L2, double* Vert_Y);


	void Setgrau(int g);
	void SetNumPto(int num);
	void SetVar(int num);
	void SetX1(double num);
	void SetX2(double num);
	void SetX3(double num);
	void SetY1(double num);
	void SetY2(double num);
	void SetY3(double num);
	void SetPesoGauss(double Peso, int& N1);
	void SetPtoInteg(double Pto, int&     N1, int&   N2);
	void SetMatriz_Jacob(double num, int N1, int N2);
	void SetMatriz_Jacob_Inv(double num, int N1, int N2);
	void SetDeterminante(double& det);
	void SetResultado(double result);

	int Getgrau();
	double* GetPesoGauss();
	double** GetPtoInteg();
	int GetNumPto();
	int GetVar();
	double GetX1();
	double GetX2();
	double GetX3();
	double GetY1();
	double GetY2();
	double GetY3();
	double GetMatriz_Jacob(int N1, int N2);
	double GetMatriz_Jacob_Inv(int N1, int N2);
	double GetDeterminante();
	
	
	void AlocPtoInteg(int&   N1, int&  N2);
	void AlocPesoGauss(int& N1);
	
	void DalocPtoInteg(int& N1);
	void DalocPesoGauss();

private:
	int grau;
	int NumPto;
	int Var;
	double**       PtoInteg;
	double*		   PesoGauss;
	double		   X1, X2, X3, Y1, Y2, Y3, Matriz_Jacob[2][2], Matriz_Jacob_Inv[2][2], Determinante;
	double		   Resultado;
	
	

};

#endif