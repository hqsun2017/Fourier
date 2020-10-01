#include "xCorr2D.h"

//循环相关
void CxCorr2D_s::CircleProcess(int nw,int nh,float* ptx,float* pty,float* ptc)
{
	int nx,ny,noff1,noff2;
	float re1,im1,*ptre1=new float[nw*nh],*ptim1=new float[nw*nh];
	float re2,im2,*ptre2=new float[nw*nh],*ptim2=new float[nw*nh];

	fft.FFT2D(nw,nh,ptx,pty,ptre1,ptim1);
	for(ny=0;ny<nh;ny++) for(nx=0;nx<nw;nx++)
	{
		noff1=ny*nw+nx;
		noff2=((nh-ny)%nh)*nw+(nw-nx)%nw;
		re1=(ptre1[noff1]+ptre1[noff2])/2;
		im1=(ptim1[noff1]-ptim1[noff2])/2;
		re2=(ptim1[noff1]+ptim1[noff2])/2;
		im2=(ptre1[noff2]-ptre1[noff1])/2;
		ptre2[noff1]=re1*re2+im1*im2;
		ptim2[noff1]=re1*im2-re2*im1;
	}
	fft.FFT2D(nw,nh,ptre2,ptim2,ptc,ptim1,false);

	delete[] ptre1; delete[] ptim1;
	delete[] ptre2; delete[] ptim2;
}

//循环相关
void CxCorr2D::CircleProcess(int nw,int nh,double* ptx,double* pty,double* ptc)
{
	int nx,ny,noff1,noff2;
	double re1,im1,*ptre1=new double[nw*nh],*ptim1=new double[nw*nh];
	double re2,im2,*ptre2=new double[nw*nh],*ptim2=new double[nw*nh];

	fft.FFT2D(nw,nh,ptx,pty,ptre1,ptim1);
	for(ny=0;ny<nh;ny++) for(nx=0;nx<nw;nx++)
	{
		noff1=ny*nw+nx;
		noff2=((nh-ny)%nh)*nw+(nw-nx)%nw;
		re1=(ptre1[noff1]+ptre1[noff2])/2;
		im1=(ptim1[noff1]-ptim1[noff2])/2;
		re2=(ptim1[noff1]+ptim1[noff2])/2;
		im2=(ptre1[noff2]-ptre1[noff1])/2;
		ptre2[noff1]=re1*re2+im1*im2;
		ptim2[noff1]=re1*im2-re2*im1;
	}
	fft.FFT2D(nw,nh,ptre2,ptim2,ptc,ptim1,false);

	delete[] ptre1; delete[] ptim1;
	delete[] ptre2; delete[] ptim2;
}
