//xPoisson2D.h : main header file

#if !defined(_XPOISSON2D_H_)
#define _XPOISSON2D_H_

#include "xFFT2D.h"
#include "xDCT_DST.h"

//Function to solve Poisson's equation in 2-d:
//d^2 P/dx^2+d^2 P/dy^2=v for xl<=x<=xh and 0<=y<=yh
//alphaL P+betaL dP/dx=gammaL(y) at x=xl
//alphaH P+betaH dP/dx=gammaH(y) at x=xh

//In y-direction,either simple Dirichlet boundary conditions: P(x,0)=P(x,yh)=0
//or simple Neumann boundary conditions: dP/dy(x,0)=dP/dy(x,yh)=0

//Matrices P and v assumed to be of extent N+2,J+1;
//Arrays gammaL,gammaH assumed to be of extent J+1.
//Now,(i,j)th elements of matrices correspond to
//x_i=xl+i*dx i=0,N+1;y_j=yl+j*dy j=0,J
//Here,dx=(xh-xl)/(N+1) dy=yh/J are grid spacing in x&y-direction.
class CxPoisson2D_s
{
public:
	CxPoisson2D_s();
	~CxPoisson2D_s() {;}

private:
	CxFFT2D_s fft2d;
	CxDCT_DST_s tran;

	int m_nWidth,m_nHeight;
	float m_tdx,m_tdy,*m_ptSource,*m_ptSolver;
	float m_txSuperior,m_txInferior,m_tySuperior;
	//Coefficients at x=xl,xh
	float m_tAlfaL,m_tBetaL,m_tAlfaH,m_tBetaH;
	float *m_ptBondL,*m_ptBondH;

private:
	void Tridiagonal(int ndim,float* pta,float* ptb,float* ptc,float* ptw,float*& ptu);

public:
	int& Width() {return m_nWidth;}
	int& Height() {return m_nHeight;}

	float& DeltaX() {return m_tdx;}
	float& DeltaY() {return m_tdy;}
	float& SupX() {return m_txSuperior;}
	float& InfX() {return m_txInferior;}
	float& SupY() {return m_tySuperior;}

	float*& Source() {return m_ptSource;}
	float*& Solver() {return m_ptSolver;}

	float& AlfaL() {return m_tAlfaL;}
	float& BetaL() {return m_tBetaL;}
	float& AlfaH() {return m_tAlfaH;}
	float& BetaH() {return m_tBetaH;}
	float*& BondL() {return m_ptBondL;}
	float*& BondH() {return m_ptBondH;}

	bool Check1();
	bool Check2();
	bool Periodic();
	bool EasyDirichlet();
	bool EasyNeumann();
	bool FastDirichlet();
	bool FastNeumann();
};

class CxPoisson2D
{
public:
	CxPoisson2D();
	~CxPoisson2D() {;}

private:
	CxFFT2D fft2d;
	CxDCT_DST tran;

	int m_nWidth,m_nHeight;
	double m_tdx,m_tdy,*m_ptSource,*m_ptSolver;
	double m_txSuperior,m_txInferior,m_tySuperior;
	//Coefficients at x=xl,xh
	double m_tAlfaL,m_tBetaL,m_tAlfaH,m_tBetaH;
	double *m_ptBondL,*m_ptBondH;

private:
	void Tridiagonal(int ndim,double* pta,double* ptb,double* ptc,double* ptw,double*& ptu);

public:
	int& Width() {return m_nWidth;}
	int& Height() {return m_nHeight;}

	double& DeltaX() {return m_tdx;}
	double& DeltaY() {return m_tdy;}
	double& SupX() {return m_txSuperior;}
	double& InfX() {return m_txInferior;}
	double& SupY() {return m_tySuperior;}

	double*& Source() {return m_ptSource;}
	double*& Solver() {return m_ptSolver;}

	double& AlfaL() {return m_tAlfaL;}
	double& BetaL() {return m_tBetaL;}
	double& AlfaH() {return m_tAlfaH;}
	double& BetaH() {return m_tBetaH;}
	double*& BondL() {return m_ptBondL;}
	double*& BondH() {return m_ptBondH;}

	bool Check1();
	bool Check2();
	bool Periodic();
	bool EasyDirichlet();
	bool EasyNeumann();
	bool FastDirichlet();
	bool FastNeumann();
};

#endif //!defined(_XPOISSON2D_H_)
