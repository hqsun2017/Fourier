//xDiffusion2D.h : main header file

#if !defined(_XDIFFUSION2D_H_)
#define _XDIFFUSION2D_H_

#include "xDCT_DST.h"

//Function to evolve diffusion equation in 2-d:
//dT/dt=D d^2 float/dx^2+D d^2 float/dy^2 for xl<=x<=xh	and 0<=y<=yh
//alphaL float+betaL dT/dx=gammaL(y) at x=xl;alphaH float+betaH dT/dx=gammaH(y) at x=xh

//In y-direction,either simple Dirichlet boundary conditions: float(x,0)=float(x,yh)=0
//or simple Neumann boundary conditions: dT/dy(x,0)=dT/dy(x,yh)=0

//Matrix float assumed to be of extent N+2,J+1;Arrays gammaL,gammaH assumed to be of extent J+1.

//Now,(i,j)th elements of matrices correspond to
//x_i=xl+i*dx i=0,N+1;y_j=yl+j*dy j=0,J
//Here,dx=(xh-xl)/(N+1) dy=yh/J are grid spacing in x&y-direction.

//Now,C=D dt/dx^2,and kappa=pi*dx/yh

//Finally,Uses Crank-Nicholson scheme.
class CxDiffusion2D_s
{
public:
	CxDiffusion2D_s();
	~CxDiffusion2D_s() {;}

private:
	CxDCT_DST_s tran;

	int m_nWidth,m_nHeight;
	float m_tD,m_tdt,m_tdx,m_tdy,*m_ptSource,*m_ptSolver;
	float m_txSuperior,m_txInferior,m_tySuperior;
	//Coefficients at x=xl,xh
	float m_tAlfaL,m_tBetaL,m_tAlfaH,m_tBetaH,*m_ptBondL,*m_ptBondH;

private:
	void Tridiagonal(int ndim,float* pta,float* ptb,float* ptc,float* ptw,float*& ptu);

public:
	int& Width() {return m_nWidth;}
	int& Height() {return m_nHeight;}

	float& CoeffD() {return m_tD;}
	float& DeltaT() {return m_tdt;}
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
	bool Dirichlet();
	bool Neumann();
};

class CxDiffusion2D
{
public:
	CxDiffusion2D();
	~CxDiffusion2D() {;}

private:
	CxDCT_DST tran;

	int m_nWidth,m_nHeight;
	double m_tD,m_tdt,m_tdx,m_tdy,*m_ptSource,*m_ptSolver;
	double m_txSuperior,m_txInferior,m_tySuperior;
	//Coefficients at x=xl,xh
	double m_tAlfaL,m_tBetaL,m_tAlfaH,m_tBetaH,*m_ptBondL,*m_ptBondH;

private:
	void Tridiagonal(int ndim,double* pta,double* ptb,double* ptc,double* ptw,double*& ptu);

public:
	int& Width() {return m_nWidth;}
	int& Height() {return m_nHeight;}

	double& CoeffD() {return m_tD;}
	double& DeltaT() {return m_tdt;}
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
	bool Dirichlet();
	bool Neumann();
};

#endif //!defined(_XDIFFUSION2D_H_)
