//xCorr1D.h : main header file

#if !defined(_XCORR1D_H_)
#define _XCORR1D_H_

#include "xFFT.h"
#include "xC2Base1D.h"

class CxCorr1D_s : public CxC2Base1D_s
{
public:
	CxCorr1D_s() {;}
	~CxCorr1D_s() {;}

private:
	CxFFT_s fft;

public:
	//循环相关
	void CircleProcess(int nSize,float* ptx,float* pty,float* ptc);

	//线性相关在纯虚类CxC2Base1D中已经定义
	//void LinearProcess(int nSizx,float* ptx,int nSizy,float* pty,float* ptc);
};

class CxCorr1D : public CxC2Base1D
{
public:
	CxCorr1D() {;}
	~CxCorr1D() {;}

private:
	CxFFT fft;

public:
	//循环相关
	void CircleProcess(int nSize,double* ptx,double* pty,double* ptc);

	//线性相关在纯虚类CxC2Base1D中已经定义
	//void LinearProcess(int nSizx,double* ptx,int nSizy,double* pty,double* ptc);
};

#endif //!defined(_XCORR1D_H_)
