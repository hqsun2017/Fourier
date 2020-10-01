//xFFT2D.h : main header file

#if !defined(_XRADIX2D_H_)
#define _XRADIX2D_H_

#include "xFFT.h"

class CxFFT2D_s
{
private:
	CxFFT_s fft;

public:
	CxFFT2D_s() {;}
	~CxFFT2D_s() {;}

public:
	void FFT2D(int nX,int nY,float* xRe,float* xIm,float* yRe,float* yIm,bool bTrans=true);
};

class CxFFT2D
{
private:
	CxFFT fft;

public:
	CxFFT2D() {;}
	~CxFFT2D() {;}

public:
	void FFT2D(int nX,int nY,double* xRe,double* xIm,double* yRe,double* yIm,bool bTrans=true);
};

#endif //!defined(_XRADIX2D_H_)
