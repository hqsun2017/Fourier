//xHilbert.h : main header file

#if !defined(_XHILBERT_H_)
#define _XHILBERT_H_

#include "xFFT.h"

class CxHilbert_s
{
public:
	CxHilbert_s() {;}
	~CxHilbert_s() {;}

private:
	CxFFT_s fft;

protected:
	//下面几种处理方法中,对称延拓的效果最好:Hilbert1
	//不作任何处理,直接利用Fourier变换对数据进行Hilbert变换
	void Hilbert0(int nSize,float* ptx,float* pty);
	//对称延拓数据序列后,利用Fourier变换对数据进行Hilbert变换
	void Hilbert1(int nSize,float* ptx,float* pty);
	//反对称延拓数据序列后,利用Fourier变换对数据进行Hilbert变换
	void Hilbert2(int nSize,float* ptx,float* pty);
	//等长补零延拓数据序列后,利用Fourier变换对数据进行Hilbert变换
	void Hilbert3(int nSize,float* ptx,float* pty);
	//加窗处理后,再利用Fourier变换对数据进行Hilbert变换
	void Hilbert4(int nSize,float* ptx,float* pty,float* ptw);

public:
	void HilbertTransform(int nSize,float* ptx,float* pty,int nOption=0,float* ptw=NULL);
};

class CxHilbert
{
public:
	CxHilbert() {;}
	~CxHilbert() {;}

private:
	CxFFT fft;

protected:
	//下面几种处理方法中,对称延拓的效果最好:Hilbert1
	//不作任何处理,直接利用Fourier变换对数据进行Hilbert变换
	void Hilbert0(int nSize,double* ptx,double* pty);
	//对称延拓数据序列后,利用Fourier变换对数据进行Hilbert变换
	void Hilbert1(int nSize,double* ptx,double* pty);
	//反对称延拓数据序列后,利用Fourier变换对数据进行Hilbert变换
	void Hilbert2(int nSize,double* ptx,double* pty);
	//等长补零延拓数据序列后,利用Fourier变换对数据进行Hilbert变换
	void Hilbert3(int nSize,double* ptx,double* pty);
	//加窗处理后,再利用Fourier变换对数据进行Hilbert变换
	void Hilbert4(int nSize,double* ptx,double* pty,double* ptw);

public:
	void HilbertTransform(int nSize,double* ptx,double* pty,int nOption=0,double* ptw=NULL);
};

#endif //!defined(_XHILBERT_H_)
