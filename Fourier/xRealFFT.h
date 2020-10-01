//xRealFFT.h : main header file

#if !defined(_XREALFFT_H_)
#define _XREALFFT_H_

#include "xFFT.h"

class CxRealFFT_s : public CxFFT_s
{
public:
	CxRealFFT_s() {;}
	~CxRealFFT_s() {;}

public:
	void RealFFT(int n,float* x1,float* x2,float* y1Re,float* y1Im,float* y2Re,float* y2Im);
};

class CxRealFFT : public CxFFT
{
public:
	CxRealFFT() {;}
	~CxRealFFT() {;}

public:
	void RealFFT(int n,double* x1,double* x2,double* y1Re,double* y1Im,double* y2Re,double* y2Im);
};

#endif //!defined(_XREALFFT_H_)
