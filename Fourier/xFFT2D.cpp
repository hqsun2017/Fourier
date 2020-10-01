#include "xFFT2D.h"

void CxFFT2D_s::FFT2D(int nX,int nY,float* xRe,float* xIm,float* yRe,float* yIm,bool bTrans)
{
	float *pxre=xRe,*pxim=xIm,*pyre=yRe,*pyim=yIm;
	for(int y=0;y<nY;y++)
	{
		fft.FFT(nX,pxre,pxim,pyre,pyim,bTrans);
		pyre+=nX; pyim+=nX;
		pxre+=nX; if(NULL!=xIm) pxim+=nX;
	}

	pxre=new float[nY]; pxim=new float[nY];
	pyre=new float[nY]; pyim=new float[nY];
	for(int x=0;x<nX;x++)
	{
		for(y=0;y<nY;y++)
		{
			pxre[y]=yRe[y*nX+x];
			pxim[y]=yIm[y*nX+x];
		}
		fft.FFT(nY,pxre,pxim,pyre,pyim,bTrans);
		for(y=0;y<nY;y++)
		{
			yRe[y*nX+x]=pyre[y];
			yIm[y*nX+x]=pyim[y];
		}
	}

	delete[] pxre; delete[] pxim;
	delete[] pyre; delete[] pyim;
}

void CxFFT2D::FFT2D(int nX,int nY,double* xRe,double* xIm,double* yRe,double* yIm,bool bTrans)
{
	double *pxre=xRe,*pxim=xIm,*pyre=yRe,*pyim=yIm;
	for(int y=0;y<nY;y++)
	{
		fft.FFT(nX,pxre,pxim,pyre,pyim,bTrans);
		pyre+=nX; pyim+=nX;
		pxre+=nX; if(NULL!=xIm) pxim+=nX;
	}

	pxre=new double[nY]; pxim=new double[nY];
	pyre=new double[nY]; pyim=new double[nY];
	for(int x=0;x<nX;x++)
	{
		for(y=0;y<nY;y++)
		{
			pxre[y]=yRe[y*nX+x];
			pxim[y]=yIm[y*nX+x];
		}
		fft.FFT(nY,pxre,pxim,pyre,pyim,bTrans);
		for(y=0;y<nY;y++)
		{
			yRe[y*nX+x]=pyre[y];
			yIm[y*nX+x]=pyim[y];
		}
	}

	delete[] pxre; delete[] pxim;
	delete[] pyre; delete[] pyim;
}
