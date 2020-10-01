//x_Corr1D.h : main header file

#if !defined(_X_CORR1D_H_)
#define _X_CORR1D_H_

#include "xConv1D.h"

class Cx_Corr1D_s : public CxConv1D_s
{
public:
	Cx_Corr1D_s() {;}
	~Cx_Corr1D_s() {;}

public:
	//循环相关
	void CircleCorr1D(int nSize,float* ptx,float* pty,float* ptc);

	//线性相关
	void LinearCorr1D(int nSizx,float* ptx,int nSizy,float* pty,float* ptc);
};

class Cx_Corr1D : public CxConv1D
{
public:
	Cx_Corr1D() {;}
	~Cx_Corr1D() {;}

public:
	//循环相关
	void CircleCorr1D(int nSize,double* ptx,double* pty,double* ptc);

	//线性相关
	void LinearCorr1D(int nSizx,double* ptx,int nSizy,double* pty,double* ptc);
};

#endif //!defined(_XCORR1D__H_)
