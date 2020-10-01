//xConv2D.h : main header file

#if !defined(_XCONV2D_H_)
#define _XCONV2D_H_

#include "xC2Base2D.h"
#include "xFFT2D.h"

class CxConv2D_s : public CxC2Base2D_s
{
public:
	CxConv2D_s() {;}
	~CxConv2D_s() {;}

private:
	CxFFT2D_s fft;

public:
	//ѭ�����
	void CircleProcess(int nw,int nh,float* ptx,float* pty,float* ptc);

	//���Ծ���ڴ�����CxC2Base2D���Ѿ�����
	//void LinearProcess(int nw1,int nh1,float* ptx,int nw2,int nh2,float* pty,float* ptc);
};

class CxConv2D : public CxC2Base2D
{
public:
	CxConv2D() {;}
	~CxConv2D() {;}

private:
	CxFFT2D fft;

public:
	//ѭ�����
	void CircleProcess(int nw,int nh,double* ptx,double* pty,double* ptc);

	//���Ծ���ڴ�����CxC2Base2D���Ѿ�����
	//void LinearProcess(int nw1,int nh1,double* ptx,int nw2,int nh2,double* pty,double* ptc);
};

#endif //!defined(_XCONV2D_H_)
