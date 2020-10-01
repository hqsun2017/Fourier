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
	//���漸�ִ�������,�Գ����ص�Ч�����:Hilbert1
	//�����κδ���,ֱ������Fourier�任�����ݽ���Hilbert�任
	void Hilbert0(int nSize,float* ptx,float* pty);
	//�Գ������������к�,����Fourier�任�����ݽ���Hilbert�任
	void Hilbert1(int nSize,float* ptx,float* pty);
	//���Գ������������к�,����Fourier�任�����ݽ���Hilbert�任
	void Hilbert2(int nSize,float* ptx,float* pty);
	//�ȳ����������������к�,����Fourier�任�����ݽ���Hilbert�任
	void Hilbert3(int nSize,float* ptx,float* pty);
	//�Ӵ������,������Fourier�任�����ݽ���Hilbert�任
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
	//���漸�ִ�������,�Գ����ص�Ч�����:Hilbert1
	//�����κδ���,ֱ������Fourier�任�����ݽ���Hilbert�任
	void Hilbert0(int nSize,double* ptx,double* pty);
	//�Գ������������к�,����Fourier�任�����ݽ���Hilbert�任
	void Hilbert1(int nSize,double* ptx,double* pty);
	//���Գ������������к�,����Fourier�任�����ݽ���Hilbert�任
	void Hilbert2(int nSize,double* ptx,double* pty);
	//�ȳ����������������к�,����Fourier�任�����ݽ���Hilbert�任
	void Hilbert3(int nSize,double* ptx,double* pty);
	//�Ӵ������,������Fourier�任�����ݽ���Hilbert�任
	void Hilbert4(int nSize,double* ptx,double* pty,double* ptw);

public:
	void HilbertTransform(int nSize,double* ptx,double* pty,int nOption=0,double* ptw=NULL);
};

#endif //!defined(_XHILBERT_H_)
