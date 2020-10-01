#include "xCorr1D.h"

void CxCorr1D_s::CircleProcess(int nSize,float* ptx,float* pty,float* ptc)
{//循环相关
	float re1,re2,im1,im2,*ptmp=new float[nSize];
	float *ptre=new float[nSize],*ptim=new float[nSize];

	fft.FFT(nSize,ptx,pty,ptre,ptim);
	int i=1,nhalf=nSize/2;
	ptre[0]=ptre[0]*ptim[0]; ptim[0]=0;
	for(;i<=nhalf;i++)
	{
		re1=(ptre[i]+ptre[nSize-i])/2;
		im1=(ptim[i]-ptim[nSize-i])/2;
		re2=(ptim[i]+ptim[nSize-i])/2;
		im2=(ptre[nSize-i]-ptre[i])/2;
		ptre[i]=re1*re2+im1*im2;
		ptim[i]=re1*im2-re2*im1;
	}
	for(;i<nSize;i++)
	{
		ptre[i]=ptre[nSize-i];
		ptim[i]=-ptim[nSize-i];
	}
	fft.FFT(nSize,ptre,ptim,ptc,ptmp,false);

	delete[] ptre; delete[] ptim; delete[] ptmp;
}

void CxCorr1D::CircleProcess(int nSize,double* ptx,double* pty,double* ptc)
{//循环相关
	double re1,re2,im1,im2,*ptmp=new double[nSize];
	double *ptre=new double[nSize],*ptim=new double[nSize];

	fft.FFT(nSize,ptx,pty,ptre,ptim);
	int i=1,nhalf=nSize/2;
	ptre[0]=ptre[0]*ptim[0]; ptim[0]=0;
	for(;i<=nhalf;i++)
	{
		re1=(ptre[i]+ptre[nSize-i])/2;
		im1=(ptim[i]-ptim[nSize-i])/2;
		re2=(ptim[i]+ptim[nSize-i])/2;
		im2=(ptre[nSize-i]-ptre[i])/2;
		ptre[i]=re1*re2+im1*im2;
		ptim[i]=re1*im2-re2*im1;
	}
	for(;i<nSize;i++)
	{
		ptre[i]=ptre[nSize-i];
		ptim[i]=-ptim[nSize-i];
	}
	fft.FFT(nSize,ptre,ptim,ptc,ptmp,false);

	delete[] ptre; delete[] ptim; delete[] ptmp;
}
