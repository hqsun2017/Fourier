//xFHT2D.h : main header file

#if !defined(_XFHT2D_H_)
#define _XFHT2D_H_

#include "xFHT.h"

class CxFHT2D_s : public CxFHT_s
{
public:
	CxFHT2D_s() {;}
	~CxFHT2D_s() {;}

public:
	bool FHT2D(int nX,int nY,float* ptHT);
	bool Correlation(int nX,int nY,float* ptx,float* pty,float* ptc);
};

class CxFHT2D : public CxFHT
{
public:
	CxFHT2D() {;}
	~CxFHT2D() {;}

public:
	bool FHT2D(int nX,int nY,double* ptHT);
	bool Correlation(int nX,int nY,double* ptx,double* pty,double* ptc);
};

#endif //!defined(_XFHT2D_H_)
