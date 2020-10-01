//xHilbert2D.h : main header file

#if !defined(_XHILBERT2D_H_)
#define _XHILBERT2D_H_

#include "xFFT2D.h"

enum HT2dMode{transform,modul,phase,orientation,RTCOS,RTSIN};

class CxHilbert2D_s : public CxFFT2D_s
{
private:

public:
	CxHilbert2D_s() {;}
	~CxHilbert2D_s() {;}

public:
	void Hilbert2D(int nWidth,int nHeight,float* ptIn,float* ptOut,int mode=transform);
	void Riesz(int nWidth,int nHeight,float* ptIn,float* ptOut,int mode=transform);
};

class CxHilbert2D : public CxFFT2D
{
private:

public:
	CxHilbert2D() {;}
	~CxHilbert2D() {;}

public:
	void Hilbert2D(int nWidth,int nHeight,double* ptIn,double* ptOut,int mode=transform);
	void Riesz(int nWidth,int nHeight,double* ptIn,double* ptOut,int mode=transform);
};

#endif //!defined(_XHILBERT2D_H_)
