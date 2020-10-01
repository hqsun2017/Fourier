//tFunc.h : main header file

#if !defined(_TFUNC_H_)
#define _TFUNC_H_

#include <math.h>
#include <stdlib.h>

void DFT_s(int N,float* pxRe,float* pxIm,float* pyRe,float* pyIm,bool bTrans=true);

void DFT2D_s(int nX,int nY,float* pxRe,float* pxIm,float* pyRe,float* pyIm,bool bTrans=true);

//HRFT:A C Program For Spectrum Magnification: When The FFT Is Not Enough,Version 1.0,March 27,2003
//by Stefan Hollos,stefan@exstrom.com;Richard Hollos,richard@exstrom.com
//Exstrom Laboratories LLC,P.O. Box 7651,Longmont,CO 80501
//Calculates the Fourier transform of nLen samples of a signal at the frequency of tOmega.
//The values of cos(n*tOmega) & sin(n*tOmega) are calculated using Chebyshev polynomials.
void HRFT_s(int nLen,float* ptIn,float tOmega,float& fRe,float& fIm);

void Rectangular_s(int N,float* ptwin);

void Triangular_s(int N,float* ptwin);

void Hanning_s(int N,float* ptwin);

void Hamming_s(int N,float* ptwin);

void Gauss_s(int N,float tpara,float* ptwin);

void Bartlett_Hanning_s(int N,float* ptwin);

void Blackman_s(int N,float* ptwin);

void Nuttall_s(int N,float* ptwin);

void Blackman_Harris_s(int N,float* ptwin);

void Blackman_Nuttall_s(int N,float* ptwin);

void Welch_s(int N,float* ptwin);

void Flat_top_s(int N,float* ptwin);

void Bohman_s(int N,float* ptwin);

//BesselI0 -- Regular Modified Cylindrical Bessel Function (Bessel I).
float BesselI0_s(float x);

void Kaiser_s(int N,float alpha,float* ptwin);

void Parzen_s(int N,float* ptwin);

void sgn_s(int nSize,float*& psgnRe,float*& psgnIm);

void SlowDCT_s(int n,float* ptx,float* pty,bool bdir=true);

void SlowDST_s(int n,float* ptx,float* pty,bool bdir=true);

void SlowDST0_s(int n,float* ptx,float* pty,bool bdir=true);

//the following definition is for the double type
void DFT(int N,double* pxRe,double* pxIm,double* pyRe,double* pyIm,bool bTrans=true);

void DFT2D(int nX,int nY,double* pxRe,double* pxIm,double* pyRe,double* pyIm,bool bTrans=true);

void HRFT(int nLen,double* ptIn,double tOmega,double& fRe,double& fIm);

void Rectangular(int N,double* ptwin);

void Triangular(int N,double* ptwin);

void Hanning(int N,double* ptwin);

void Hamming(int N,double* ptwin);

void Gauss(int N,double tpara,double* ptwin);

void Bartlett_Hanning(int N,double* ptwin);

void Blackman(int N,double* ptwin);

void Nuttall(int N,double* ptwin);

void Blackman_Harris(int N,double* ptwin);

void Blackman_Nuttall(int N,double* ptwin);

void Welch(int N,double* ptwin);

void Flat_top(int N,double* ptwin);

void Bohman(int N,double* ptwin);

double BesselI0(double x);

void Kaiser(int N,double alpha,double* ptwin);

void Parzen(int N,double* ptwin);

void sgn(int nSize,double*& psgnRe,double*& psgnIm);

void SlowDCT(int n,double* ptx,double* pty,bool bdir=true);

void SlowDST(int n,double* ptx,double* pty,bool bdir=true);

void SlowDST0(int n,double* ptx,double* pty,bool bdir=true);

#endif //!defined(_TFUNC_H_)
