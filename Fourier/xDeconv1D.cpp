#include "xDeconv1D.h"

//循环解卷积
int CxDeconv1D_s::Deconvolution(int nSize,float* ptIn,float* ptOut,float* ptSys)
{
	float sqr,re1,re2,im1,im2,*ptre=new float[nSize],*ptim=new float[nSize],*ptmp=new float[nSize];

	fft.FFT(nSize,ptIn,ptOut,ptre,ptim);
	int nzeros=0,i=1,nhalf=nSize/2;
	sqr=ptre[0]*ptre[0];
	if(PRECISION>sqr) {nzeros++; sqr=(float)PRECISION;}
	ptre[0]=ptre[0]*ptim[0]/sqr; ptim[0]=0;
	for(;i<=nhalf;i++)
	{
		re1=(ptre[i]+ptre[nSize-i])/2;
		im1=(ptim[i]-ptim[nSize-i])/2;
		re2=(ptim[i]+ptim[nSize-i])/2;
		im2=(ptre[nSize-i]-ptre[i])/2;
		sqr=re1*re1+im1*im1;
		if(PRECISION>sqr) {nzeros++; sqr=(float)PRECISION;}
		ptre[i]=(re1*re2+im1*im2)/sqr;
		ptim[i]=(re1*im2-re2*im1)/sqr;
	}
	for(;i<nSize;i++)
	{
		ptre[i]=ptre[nSize-i];
		ptim[i]=-ptim[nSize-i];
	}
	fft.FFT(nSize,ptre,ptim,ptSys,ptmp,false);

	delete[] ptre; delete[] ptim; delete[] ptmp;

	return nzeros;
}

//循环解卷积
int CxDeconv1D::Deconvolution(int nSize,double* ptIn,double* ptOut,double* ptSys)
{
	double sqr,re1,re2,im1,im2,*ptre=new double[nSize],*ptim=new double[nSize],*ptmp=new double[nSize];

	fft.FFT(nSize,ptIn,ptOut,ptre,ptim);
	int nzeros=0,i=1,nhalf=nSize/2;
	sqr=ptre[0]*ptre[0];
	if(PRECISION>sqr) {nzeros++; sqr=PRECISION;}
	ptre[0]=ptre[0]*ptim[0]/sqr; ptim[0]=0;
	for(;i<=nhalf;i++)
	{
		re1=(ptre[i]+ptre[nSize-i])/2;
		im1=(ptim[i]-ptim[nSize-i])/2;
		re2=(ptim[i]+ptim[nSize-i])/2;
		im2=(ptre[nSize-i]-ptre[i])/2;
		sqr=re1*re1+im1*im1;
		if(PRECISION>sqr) {nzeros++; sqr=PRECISION;}
		ptre[i]=(re1*re2+im1*im2)/sqr;
		ptim[i]=(re1*im2-re2*im1)/sqr;
	}
	for(;i<nSize;i++)
	{
		ptre[i]=ptre[nSize-i];
		ptim[i]=-ptim[nSize-i];
	}
	fft.FFT(nSize,ptre,ptim,ptSys,ptmp,false);

	delete[] ptre; delete[] ptim; delete[] ptmp;

	return nzeros;
}
