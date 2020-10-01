//xFFT.h : main header file

#if !defined(_PRIMEFACTORFFT_H_)
#define _PRIMEFACTORFFT_H_

#if !defined(NULL)
#define NULL 0
#endif

#include <math.h>

#define PRIMEFACTOR 37
//PRIMEFACTOR_2=(PRIMEFACTOR+1)/2
#define PRIMEFACTOR_2 19
#define PRIMECOUNT 20

class CxFFT_s
{
private:
	float tPI,t2PI,tIRT2,tIRT3;//t2PI=2*PI,tIRT2=1/sqrt(2),tIRT3=1/sqrt(3)

	//Variables for PrimeFFT
	bool bAlloc;
	int m_nOldSize,m_nFactor,nNewFactorSize,nHalfFactorSize;
	int m_sofarRadix[PRIMECOUNT],m_actualRadix[PRIMECOUNT],m_remainRadix[PRIMECOUNT];

	float c3_1;
	float c5_1,c5_2,c5_3,c5_4,c5_5;
	float c7_1,c7_2,c7_3,c7_4,c7_5,c7_6;
	float c9_2,c9_3,c9_4,c9_5,c9_6,c9_7,c9_8,c9_9;
	float c11_1,c11_2,c11_3,c11_4,c11_5,c11_6,c11_7,c11_8,c11_9,c11_10;
	float c13_1,c13_2,c13_3,c13_4,c13_5,c13_6,c13_7,c13_8,c13_9,c13_10,c13_11,c13_12;
	float c16_2,c16_3,c16_4,c16_5;

	int groupOffset,dataOffset,blockOffset;
	int adr,groupNo,dataNo,blockNo,twNo;
	float omega,tw_re,tw_im;
	float *pftwRe,*pftgRe,*pfzRe,*pfvRe,*pfwRe;
	float *pftwIm,*pftgIm,*pfzIm,*pfvIm,*pfwIm;
	float twiddleRe[PRIMEFACTOR],trigRe[PRIMEFACTOR],zRe[PRIMEFACTOR];
	float twiddleIm[PRIMEFACTOR],trigIm[PRIMEFACTOR],zIm[PRIMEFACTOR];
	float vRe[PRIMEFACTOR_2],wRe[PRIMEFACTOR_2];
	float vIm[PRIMEFACTOR_2],wIm[PRIMEFACTOR_2];

	//Temporary variables 
	float ttmp;
	float t1_re,t1_im,t2_re,t2_im,t3_re,t3_im,t4_re,t4_im,t5_re,t5_im,t6_re,t6_im;
	float t7_re,t7_im,t8_re,t8_im,t9_re,t9_im,t10_re,t10_im,t11_re,t11_im,t12_re,t12_im;
	float t13_re,t13_im,t14_re,t14_im,t15_re,t15_im,t16_re,t16_im,t17_re,t17_im;
	float t18_re,t18_im,t19_re,t19_im,t20_re,t20_im,t21_re,t21_im,t22_re,t22_im;

	float m1_re,m1_im,m2_re,m2_im,m3_re,m3_im,m4_re,m4_im,m5_re,m5_im,m6_re,m6_im;
	float m7_re,m7_im,m8_re,m8_im,m9_re,m9_im,m10_re,m10_im,m11_re,m11_im,m12_re,m12_im;

public:
	CxFFT_s();
	~CxFFT_s();

protected:
	//Member functions
	void ReleaseMem();
	void AllocateMem();
	void Factorize(int n,int& nFact,int* fact);
	void PrimeSetup(int nPoints);
	void Permute(int nPoint,float* xRe,float* xIm,float* yRe,float* yIm,bool bTrans=true);
	void InitTrig(int radix);
	void Radix_2(float* aRe,float* aIm);
	void Radix_3(float* aRe,float* aIm);
	void Radix_4(float* aRe,float* aIm);
	void Radix_5(float* aRe,float* aIm);
	void Radix_7(float* aRe,float* aIm);
	void Radix_8(float* aRe,float* aIm);
	void Radix_9(float* aRe,float* aIm);
	void Radix_10(float* aRe,float* aIm);
	void Radix_11(float* aRe,float* aIm);
	void Radix_13(float* aRe,float* aIm);
	void Radix_16(float* aRe,float* aIm);
	void Radix_Other(int radix);
	void TwiddleFFT(int sofarRadix,int radix,int remainRadix,float* yRe,float* yIm);

public:
	void FFT(int n,float* xRe,float* xIm,float* yRe,float* yIm,bool bTrans=true);
};

class CxFFT
{
private:
	double tPI,t2PI,tIRT2,tIRT3;//t2PI=2*PI,tIRT2=1/sqrt(2),tIRT3=1/sqrt(3)

	//Variables for PrimeFFT
	bool bAlloc;
	int m_nOldSize,m_nFactor,nNewFactorSize,nHalfFactorSize;
	int m_sofarRadix[PRIMECOUNT],m_actualRadix[PRIMECOUNT],m_remainRadix[PRIMECOUNT];

	double c3_1;
	double c5_1,c5_2,c5_3,c5_4,c5_5;
	double c7_1,c7_2,c7_3,c7_4,c7_5,c7_6;
	double c9_2,c9_3,c9_4,c9_5,c9_6,c9_7,c9_8,c9_9;
	double c11_1,c11_2,c11_3,c11_4,c11_5,c11_6,c11_7,c11_8,c11_9,c11_10;
	double c13_1,c13_2,c13_3,c13_4,c13_5,c13_6,c13_7,c13_8,c13_9,c13_10,c13_11,c13_12;
	double c16_2,c16_3,c16_4,c16_5;

	int groupOffset,dataOffset,blockOffset;
	int adr,groupNo,dataNo,blockNo,twNo;
	double omega,tw_re,tw_im;
	double *pftwRe,*pftgRe,*pfzRe,*pfvRe,*pfwRe;
	double *pftwIm,*pftgIm,*pfzIm,*pfvIm,*pfwIm;
	double twiddleRe[PRIMEFACTOR],trigRe[PRIMEFACTOR],zRe[PRIMEFACTOR];
	double twiddleIm[PRIMEFACTOR],trigIm[PRIMEFACTOR],zIm[PRIMEFACTOR];
	double vRe[PRIMEFACTOR_2],wRe[PRIMEFACTOR_2];
	double vIm[PRIMEFACTOR_2],wIm[PRIMEFACTOR_2];

	//Temporary variables 
	double ttmp;
	double t1_re,t1_im,t2_re,t2_im,t3_re,t3_im,t4_re,t4_im,t5_re,t5_im,t6_re,t6_im;
	double t7_re,t7_im,t8_re,t8_im,t9_re,t9_im,t10_re,t10_im,t11_re,t11_im,t12_re,t12_im;
	double t13_re,t13_im,t14_re,t14_im,t15_re,t15_im,t16_re,t16_im,t17_re,t17_im;
	double t18_re,t18_im,t19_re,t19_im,t20_re,t20_im,t21_re,t21_im,t22_re,t22_im;

	double m1_re,m1_im,m2_re,m2_im,m3_re,m3_im,m4_re,m4_im,m5_re,m5_im,m6_re,m6_im;
	double m7_re,m7_im,m8_re,m8_im,m9_re,m9_im,m10_re,m10_im,m11_re,m11_im,m12_re,m12_im;

public:
	CxFFT();
	~CxFFT();

protected:
	//Member functions
	void ReleaseMem();
	void AllocateMem();
	void Factorize(int n,int& nFact,int* fact);
	void PrimeSetup(int nPoints);
	void Permute(int nPoint,double* xRe,double* xIm,double* yRe,double* yIm,bool bTrans=true);
	void InitTrig(int radix);
	void Radix_2(double* aRe,double* aIm);
	void Radix_3(double* aRe,double* aIm);
	void Radix_4(double* aRe,double* aIm);
	void Radix_5(double* aRe,double* aIm);
	void Radix_7(double* aRe,double* aIm);
	void Radix_8(double* aRe,double* aIm);
	void Radix_9(double* aRe,double* aIm);
	void Radix_10(double* aRe,double* aIm);
	void Radix_11(double* aRe,double* aIm);
	void Radix_13(double* aRe,double* aIm);
	void Radix_16(double* aRe,double* aIm);
	void Radix_Other(int radix);
	void TwiddleFFT(int sofarRadix,int radix,int remainRadix,double* yRe,double* yIm);

public:
	void FFT(int n,double* xRe,double* xIm,double* yRe,double* yIm,bool bTrans=true);
};

#endif //!defined(_PRIMEFACTORFFT_H_)
