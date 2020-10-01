//xConv1D.h : main header file

#if !defined(_XCONV1D_H_)
#define _XCONV1D_H_

#include "xFFT.h"
#include "xC2Base1D.h"

class CxConv1D_s : public CxC2Base1D_s
{
public:
	CxConv1D_s() {;}
	~CxConv1D_s() {;}

private:
	CxFFT_s fft;

public:
	//循环卷积
	void CircleProcess(int nSize,float* ptx,float* pty,float* ptc);

	//线性卷积在纯虚类CxC2Base1D_s中已经定义
	//void LinearProcess(int nSizx,float* ptx,int nSizy,float* pty,float* ptc);

	//重叠相加法计算线性卷积
	void OverlapConv1D(int nSizx,float* ptx,int nSizy,float* pty,float* ptc);
};

class CxConv1D : public CxC2Base1D
{
public:
	CxConv1D() {;}
	~CxConv1D() {;}

private:
	CxFFT fft;

public:
	//循环卷积
	void CircleProcess(int nSize,double* ptx,double* pty,double* ptc);

	//线性卷积在纯虚类CxC2Base1D中已经定义
	//void LinearProcess(int nSizx,double* ptx,int nSizy,double* pty,double* ptc);

	//重叠相加法计算线性卷积
	void OverlapConv1D(int nSizx,double* ptx,int nSizy,double* pty,double* ptc);
};

#endif //!defined(_XCONV1D_H_)
