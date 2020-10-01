#include "xPhase2D.h"

void CxPhase2D_s::CircleProcess(int nw,int nh,float* ptx,float* pty,float* ptc)
{//循环相关
	int nx,ny;
	float re1,im1,*ptre1=new float[nw*nh],*ptim1=new float[nw*nh];
	float re2,im2,*ptre2=new float[nw*nh],*ptim2=new float[nw*nh];
	float re,im,mod;

	fft.FFT2D(nw,nh,ptx,pty,ptre1,ptim1);
	for(ny=0;ny<nh;ny++) for(nx=0;nx<nw;nx++)
	{
		re1=(ptre1[ny*nw+nx]+ptre1[((nh-ny)%nh)*nw+(nw-nx)%nw])/2;
		im1=(ptim1[ny*nw+nx]-ptim1[((nh-ny)%nh)*nw+(nw-nx)%nw])/2;
		re2=(ptim1[ny*nw+nx]+ptim1[((nh-ny)%nh)*nw+(nw-nx)%nw])/2;
		im2=(ptre1[((nh-ny)%nh)*nw+(nw-nx)%nw]-ptre1[ny*nw+nx])/2;
		re=re1*re2+im1*im2; im=re1*im2-re2*im1; mod=(float)sqrt(re*re+im*im);
		//if(0==mod) ptre2[ny*nw+nx]=ptim2[ny*nw+nx]=0;
		if(PRECISION>mod) ptre2[ny*nw+nx]=ptim2[ny*nw+nx]=0;
		else
		{
			ptre2[ny*nw+nx]=re/mod;
			ptim2[ny*nw+nx]=im/mod;
		}
	}
	fft.FFT2D(nw,nh,ptre2,ptim2,ptc,ptim1,false);

	delete[] ptre1; delete[] ptim1;
	delete[] ptre2; delete[] ptim2;
}

void CxPhase2D::CircleProcess(int nw,int nh,double* ptx,double* pty,double* ptc)
{//循环相关
	int nx,ny;
	double re1,im1,*ptre1=new double[nw*nh],*ptim1=new double[nw*nh];
	double re2,im2,*ptre2=new double[nw*nh],*ptim2=new double[nw*nh];
	double re,im,mod;

	fft.FFT2D(nw,nh,ptx,pty,ptre1,ptim1);
	for(ny=0;ny<nh;ny++) for(nx=0;nx<nw;nx++)
	{
		re1=(ptre1[ny*nw+nx]+ptre1[((nh-ny)%nh)*nw+(nw-nx)%nw])/2;
		im1=(ptim1[ny*nw+nx]-ptim1[((nh-ny)%nh)*nw+(nw-nx)%nw])/2;
		re2=(ptim1[ny*nw+nx]+ptim1[((nh-ny)%nh)*nw+(nw-nx)%nw])/2;
		im2=(ptre1[((nh-ny)%nh)*nw+(nw-nx)%nw]-ptre1[ny*nw+nx])/2;
		re=re1*re2+im1*im2; im=re1*im2-re2*im1; mod=sqrt(re*re+im*im);
		//if(0==mod) ptre2[ny*nw+nx]=ptim2[ny*nw+nx]=0;
		if(PRECISION>mod) ptre2[ny*nw+nx]=ptim2[ny*nw+nx]=0;
		else
		{
			ptre2[ny*nw+nx]=re/mod;
			ptim2[ny*nw+nx]=im/mod;
		}
	}
	fft.FFT2D(nw,nh,ptre2,ptim2,ptc,ptim1,false);

	delete[] ptre1; delete[] ptim1;
	delete[] ptre2; delete[] ptim2;
}
