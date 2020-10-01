//xFHT.h : main header file

#if !defined(_XFHT_H_)
#define _XFHT_H_

#include <math.h>

//Note:The original source code is downlodaded from
//http://www.geocities.com/ResearchTriangle/8869/fft_summary.html

#define TRIG 22

class CxFHT_s
{
private:
	float tSQRT2;//tSQRT2=sqrt(2)
	float tCosTab[TRIG],tSinTab[TRIG];

public:
	CxFHT_s();
	~CxFHT_s() {;}

protected:
	bool powerof2(int n);

public:
	bool FHT(int nlen,float* ptHT);
	bool Spectrum(int nlen,float* ptHT);
	bool Convolution(int nlen,float* ptx,float* pty,float* ptc);
	bool Correlation(int nlen,float* ptx,float* pty,float* ptc);
};

class CxFHT
{
private:
	double tSQRT2;//tSQRT2=sqrt(2)
	double tCosTab[TRIG],tSinTab[TRIG];

public:
	CxFHT();
	~CxFHT() {;}

protected:
	bool powerof2(int n);

public:
	bool FHT(int nlen,double* ptHT);
	bool Spectrum(int nlen,double* ptHT);
	bool Convolution(int nlen,double* ptx,double* pty,double* ptc);
	bool Correlation(int nlen,double* ptx,double* pty,double* ptc);
};

#endif //!defined(_XFHT_H_)
