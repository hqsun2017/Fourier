//xDeconv1D.h : main header file

#if !defined(_XDECONV1D_H_)
#define _XDECONV1D_H_

#ifndef PRECISION
#define PRECISION 1e-20
#endif

#include "xFFT.h"

class CxDeconv1D_s
{
public:
	CxDeconv1D_s() {;}
	~CxDeconv1D_s() {;}

private:
	CxFFT_s fft;

public:
	//循环解卷积
	int Deconvolution(int nSize,float* ptIn,float* ptOut,float* ptSys);
};

class CxDeconv1D
{
public:
	CxDeconv1D() {;}
	~CxDeconv1D() {;}

private:
	CxFFT fft;

public:
	//循环解卷积
	int Deconvolution(int nSize,double* ptIn,double* ptOut,double* ptSys);
};

#endif //!defined(_XDECONV1D_H_)
