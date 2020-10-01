//xCorr2D.h : main header file

#if !defined(_XCORR2D_H_)
#define _XCORR2D_H_

#include "xFFT2D.h"
#include "xC2Base2D.h"

class CxCorr2D_s : public CxC2Base2D_s
{
public:
	CxCorr2D_s() {;}
	~CxCorr2D_s() {;}

private:
	CxFFT2D_s fft;

public:
	//ѭ�����
	void CircleProcess(int nw,int nh,float* ptx,float* pty,float* ptc);

	//��������ڴ�����CxC2Base2D���Ѿ�����
	//void LinearProcess(int nw1,int nh1,float* ptx,int nw2,int nh2,float* pty,float* ptc);
};

class CxCorr2D : public CxC2Base2D
{
public:
	CxCorr2D() {;}
	~CxCorr2D() {;}

private:
	CxFFT2D fft;

public:
	//ѭ�����
	void CircleProcess(int nw,int nh,double* ptx,double* pty,double* ptc);

	//��������ڴ�����CxC2Base2D���Ѿ�����
	//void LinearProcess(int nw1,int nh1,double* ptx,int nw2,int nh2,double* pty,double* ptc);
};

#endif //!defined(_XCORR2D_H_)
