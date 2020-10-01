//xDeconv2D.h : main header file

#if !defined(_XDECONV2D_H_)
#define _XDECONV2D_H_

#ifndef PRECISION
#define PRECISION 1e-20
#endif

#include "xFFT2D.h"

class CxDeconv2D_s
{
public:
	CxDeconv2D_s() {;}
	~CxDeconv2D_s() {;}

private:
	CxFFT2D_s fft;

public:
	//循环解卷积
	int Deconvolution(int nw,int nh,float* ptIn,float* ptOut,float* ptSys);
};

class CxDeconv2D
{
public:
	CxDeconv2D() {;}
	~CxDeconv2D() {;}

private:
	CxFFT2D fft;

public:
	//循环解卷积
	int Deconvolution(int nw,int nh,double* ptIn,double* ptOut,double* ptSys);
};

#endif //!defined(_XDECONV2D_H_)
