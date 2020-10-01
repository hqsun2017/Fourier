#include "xDeconv2D.h"

int CxDeconv2D_s::Deconvolution(int nw,int nh,float* ptIn,float* ptOut,float* ptSys)
{//循环解卷积
	int nzeros=0,nx,ny,noff1,noff2;
	float re1,im1,*ptre1=new float[nw*nh],*ptim1=new float[nw*nh];
	float re2,im2,*ptre2=new float[nw*nh],*ptim2=new float[nw*nh];
	float sqr;

	fft.FFT2D(nw,nh,ptIn,ptOut,ptre1,ptim1);
	for(ny=0;ny<nh;ny++) for(nx=0;nx<nw;nx++)
	{
		noff1=ny*nw+nx;
		noff2=((nh-ny)%nh)*nw+(nw-nx)%nw;
		re1=(ptre1[noff1]+ptre1[noff2])/2;
		im1=(ptim1[noff1]-ptim1[noff2])/2;
		re2=(ptim1[noff1]+ptim1[noff2])/2;
		im2=(ptre1[noff2]-ptre1[noff1])/2;
		sqr=re1*re1+im1*im1;
		if(PRECISION>sqr) {nzeros++; sqr=(float)PRECISION;}
		ptre2[noff1]=(re1*re2+im1*im2)/sqr;
		ptim2[noff1]=(re1*im2-re2*im1)/sqr;
	}
	fft.FFT2D(nw,nh,ptre2,ptim2,ptSys,ptim1,false);

	delete[] ptre1; delete[] ptim1;
	delete[] ptre2; delete[] ptim2;

	return nzeros;
}

int CxDeconv2D::Deconvolution(int nw,int nh,double* ptIn,double* ptOut,double* ptSys)
{//循环解卷积
	int nzeros=0,nx,ny,noff1,noff2;
	double re1,im1,*ptre1=new double[nw*nh],*ptim1=new double[nw*nh];
	double re2,im2,*ptre2=new double[nw*nh],*ptim2=new double[nw*nh];
	double sqr;

	fft.FFT2D(nw,nh,ptIn,ptOut,ptre1,ptim1);
	for(ny=0;ny<nh;ny++) for(nx=0;nx<nw;nx++)
	{
		noff1=ny*nw+nx;
		noff2=((nh-ny)%nh)*nw+(nw-nx)%nw;
		re1=(ptre1[noff1]+ptre1[noff2])/2;
		im1=(ptim1[noff1]-ptim1[noff2])/2;
		re2=(ptim1[noff1]+ptim1[noff2])/2;
		im2=(ptre1[noff2]-ptre1[noff1])/2;
		sqr=re1*re1+im1*im1;
		if(PRECISION>sqr) {nzeros++; sqr=PRECISION;}
		ptre2[noff1]=(re1*re2+im1*im2)/sqr;
		ptim2[noff1]=(re1*im2-re2*im1)/sqr;
	}
	fft.FFT2D(nw,nh,ptre2,ptim2,ptSys,ptim1,false);

	delete[] ptre1; delete[] ptim1;
	delete[] ptre2; delete[] ptim2;

	return nzeros;
}
