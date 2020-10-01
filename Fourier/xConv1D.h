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
	//ѭ�����
	void CircleProcess(int nSize,float* ptx,float* pty,float* ptc);

	//���Ծ���ڴ�����CxC2Base1D_s���Ѿ�����
	//void LinearProcess(int nSizx,float* ptx,int nSizy,float* pty,float* ptc);

	//�ص���ӷ��������Ծ��
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
	//ѭ�����
	void CircleProcess(int nSize,double* ptx,double* pty,double* ptc);

	//���Ծ���ڴ�����CxC2Base1D���Ѿ�����
	//void LinearProcess(int nSizx,double* ptx,int nSizy,double* pty,double* ptc);

	//�ص���ӷ��������Ծ��
	void OverlapConv1D(int nSizx,double* ptx,int nSizy,double* pty,double* ptc);
};

#endif //!defined(_XCONV1D_H_)
