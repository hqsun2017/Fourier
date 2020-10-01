//xDCT_DST.h : main header file

#if !defined(_XDCT_DST_H_)
#define _XDCT_DST_H_

#include "xFFT.h"

class CxDCT_DST_s : public CxFFT_s
{
public:
	CxDCT_DST_s() {;}
	~CxDCT_DST_s() {;}

public:
	void EasyDCT(int n,float* ptx,float* pty,bool bdir=true);
	void EasyDST(int n,float* ptx,float* pty,bool bdir=true);
	void EasyDST0(int n,float* ptx,float* pty,bool bdir=true);

	void FastDCT(int n,float* ptx1,float* ptx2,float* pty1,float* pty2,bool bdir=true);
	void FastDST(int n,float* ptx1,float* ptx2,float* pty1,float* pty2,bool bdir=true);
	void FastDST0(int n,float* ptx1,float* ptx2,float* pty1,float* pty2,bool bdir=true);
};

class CxDCT_DST : public CxFFT
{
public:
	CxDCT_DST() {;}
	~CxDCT_DST() {;}

public:
	void EasyDCT(int n,double* ptx,double* pty,bool bdir=true);
	void EasyDST(int n,double* ptx,double* pty,bool bdir=true);
	void EasyDST0(int n,double* ptx,double* pty,bool bdir=true);

	void FastDCT(int n,double* ptx1,double* ptx2,double* pty1,double* pty2,bool bdir=true);
	void FastDST(int n,double* ptx1,double* ptx2,double* pty1,double* pty2,bool bdir=true);
	void FastDST0(int n,double* ptx1,double* ptx2,double* pty1,double* pty2,bool bdir=true);
};

#endif //!defined(_XDCT_DST_H_)
