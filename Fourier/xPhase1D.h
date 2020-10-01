//xPhase1D.h : main header file

#if !defined(_XPHASE1D_H_)
#define _XPHASE1D_H_

#include "xFFT.h"
#include "xC2Base1D.h"

#ifndef PRECISION
#define PRECISION 1e-20
#endif

class CxPhase1D_s : public CxC2Base1D_s
{
public:
	CxPhase1D_s() {;}
	~CxPhase1D_s() {;}

private:
	CxFFT_s fft;

public:
	//循环相关
	void CircleProcess(int nSize,float* ptx,float* pty,float* ptc);

	//线性卷积在纯虚类CxC2Base1D中已经定义
	//void LinearProcess(int nSizx,float* ptx,int nSizy,float* pty,float* ptc);
};

class CxPhase1D : public CxC2Base1D
{
public:
	CxPhase1D() {;}
	~CxPhase1D() {;}

private:
	CxFFT fft;

public:
	//循环相关
	void CircleProcess(int nSize,double* ptx,double* pty,double* ptc);

	//线性卷积在纯虚类CxC2Base1D中已经定义
	//void LinearProcess(int nSizx,double* ptx,int nSizy,double* pty,double* ptc);
};

#endif //!defined(_XPHASE1D_H_)
