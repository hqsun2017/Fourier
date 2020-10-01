#include "xPhase1D.h"

/*void CxPhase1D_s::CircleProcess(int nSize,float* ptx,float* pty,float* ptc)
{//循环相关
	float re1,im1,*ptre1=new float[nSize],*ptim1=new float[nSize];
	float re2,im2,*ptre2=new float[nSize],*ptim2=new float[nSize];
	float mod,re,im;

	fft.FFT(nSize,ptx,pty,ptre1,ptim1);
	for(int i=0;i<nSize;i++)
	{
		re1=(ptre1[i]+ptre1[(nSize-i)%nSize])/2;
		im1=(ptim1[i]-ptim1[(nSize-i)%nSize])/2;
		re2=(ptim1[i]+ptim1[(nSize-i)%nSize])/2;
		im2=(ptre1[(nSize-i)%nSize]-ptre1[i])/2;
		re=re1*re2+im1*im2; im=re1*im2-re2*im1;
		mod=sqrt(re*re+im*im);

		if(PRECISION>mod) ptre2[i]=ptim2[i]=0;
		else
		{
			ptre2[i]=re/mod;
			ptim2[i]=im/mod;
		}
	}
	fft.FFT(nSize,ptre2,ptim2,ptc,ptim1,false);

	delete[] ptre1; delete[] ptim1;
	delete[] ptre2; delete[] ptim2;
}*/

void CxPhase1D_s::CircleProcess(int nSize,float* ptx,float* pty,float* ptc)
{//循环相关
	float mod,re1,re2,im1,im2,*ptmp=new float[nSize];
	float re,im,*ptre=new float[nSize],*ptim=new float[nSize];

	fft.FFT(nSize,ptx,pty,ptre,ptim);
	int i=1,nhalf=nSize/2;
	ptre[0]=ptre[0]*ptim[0]; ptim[0]=0;
	if(PRECISION>fabs(ptre[0])) ptre[0]=0;
	else ptre[0]=0<ptre[0] ? 1.0f : -1.0f;
	for(;i<=nhalf;i++)
	{
		re1=(ptre[i]+ptre[nSize-i])/2;
		im1=(ptim[i]-ptim[nSize-i])/2;
		re2=(ptim[i]+ptim[nSize-i])/2;
		im2=(ptre[nSize-i]-ptre[i])/2;
		re=re1*re2+im1*im2; im=re1*im2-re2*im1;
		mod=(float)sqrt(re*re+im*im);
		if(PRECISION>mod) ptre[i]=ptim[i]=0;
		else
		{
			ptre[i]=re/mod;
			ptim[i]=im/mod;
		}
	}
	for(;i<nSize;i++)
	{
		ptre[i]=ptre[nSize-i];
		ptim[i]=-ptim[nSize-i];
	}
	fft.FFT(nSize,ptre,ptim,ptc,ptmp,false);

	delete[] ptre; delete[] ptim; delete[] ptmp;
}

void CxPhase1D::CircleProcess(int nSize,double* ptx,double* pty,double* ptc)
{//循环相关
	double mod,re1,re2,im1,im2,*ptmp=new double[nSize];
	double re,im,*ptre=new double[nSize],*ptim=new double[nSize];

	fft.FFT(nSize,ptx,pty,ptre,ptim);
	int i=1,nhalf=nSize/2;
	ptre[0]=ptre[0]*ptim[0]; ptim[0]=0;
	if(PRECISION>fabs(ptre[0])) ptre[0]=0;
	else ptre[0]=0<ptre[0] ? 1 : -1;
	for(;i<=nhalf;i++)
	{
		re1=(ptre[i]+ptre[nSize-i])/2;
		im1=(ptim[i]-ptim[nSize-i])/2;
		re2=(ptim[i]+ptim[nSize-i])/2;
		im2=(ptre[nSize-i]-ptre[i])/2;
		re=re1*re2+im1*im2; im=re1*im2-re2*im1;
		mod=sqrt(re*re+im*im);
		if(PRECISION>mod) ptre[i]=ptim[i]=0;
		else
		{
			ptre[i]=re/mod;
			ptim[i]=im/mod;
		}
	}
	for(;i<nSize;i++)
	{
		ptre[i]=ptre[nSize-i];
		ptim[i]=-ptim[nSize-i];
	}
	fft.FFT(nSize,ptre,ptim,ptc,ptmp,false);

	delete[] ptre; delete[] ptim; delete[] ptmp;
}
