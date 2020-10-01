#include "xHilbert2D.h"

void CxHilbert2D_s::Hilbert2D(int nWidth,int nHeight,float* ptIn,float* ptOut,int mode)
{
	int nx,ny,nsize=nWidth*nHeight;

	//////////////////////////////////////////////////////////////////////
	//calculate directly the 2d transform
	/*float *pxre=new float[nsize],*pxim=new float[nsize],*pyim=new float[nsize];
	FFT2D(nWidth,nHeight,ptIn,NULL,pxre,pxim);
	for(ny=0;ny<nHeight;ny++) pxim[ny*nWidth]=pxre[ny*nWidth]=0;
	for(ny=0;ny<nHeight;ny++) for(nx=1;nx<nWidth;nx++)
	{
		int noff=ny*nWidth+nx;
		if(nx>nWidth/2) pxim[noff]*=-1;
		else pxre[noff]*=-1;
	}
	FFT2D(nWidth,nHeight,pxim,pxre,ptOut,pyim,false);
	delete[] pxre; delete[] pxim; delete[] pyim;*/

	//calculate the 2d transform via 1d Hilbert transform
	CxFFT_s fft;
	float *ptFTRe=new float[nWidth],*ptFTIm=new float[nWidth],*ptHTRe=new float[nWidth];
	float *ptx=ptIn,*pty=ptOut;
	for(ny=0;ny<nHeight;ny++)
	{
		fft.FFT(nWidth,ptx,NULL,ptFTRe,ptFTIm);
		if(1==nWidth%2)
		{
			for(nx=1;nx<=(nWidth-1)/2;nx++) {ptFTRe[nx]*=2; ptFTIm[nx]*=2;}
			for(nx=(nWidth+1)/2;nx<nWidth;nx++) ptFTRe[nx]=ptFTIm[nx]=0;
		}
		else
		{
			for(nx=1;nx<=nWidth/2;nx++) {ptFTRe[nx]*=2; ptFTIm[nx]*=2;}
			for(nx=nWidth/2+1;nx<nWidth;nx++) ptFTRe[nx]=ptFTIm[nx]=0;
		}
		fft.FFT(nWidth,ptFTRe,ptFTIm,ptHTRe,pty,false);
		ptx+=nWidth; pty+=nWidth;
	}
	delete[] ptFTRe; delete[] ptFTIm; delete[] ptHTRe;
	//////////////////////////////////////////////////////////////////////

	for(nx=0;nx<nsize;nx++)
	{
		if(modul==mode) ptOut[nx]=(float)sqrt(ptIn[nx]*ptIn[nx]+ptOut[nx]*ptOut[nx]);
		else if(phase==mode) ptOut[nx]=(float)atan2(ptOut[nx],ptIn[nx]);
	}
}

void CxHilbert2D_s::Riesz(int nWidth,int nHeight,float* ptIn,float* ptOut,int mode)
{
	int nx,ny,nxp,nyp,noff,nsize=nWidth*nHeight;
	float tmp,*pfftre=new float[nsize],*pfftim=new float[nsize];
	float tcos,*pcosre=new float[nsize],*pcosim=new float[nsize];
	float tsin,*psinre=new float[nsize],*psinim=new float[nsize];

	FFT2D(nWidth,nHeight,ptIn,NULL,pfftre,pfftim);

	for(ny=0;ny<nHeight;ny++) for(nx=0;nx<nWidth;nx++)
	{
		//nxp=(nx+nWidth/2)%nWidth-nWidth/2;//nxp=(nx<nWidth/2) ? nx : nx-nWidth;
		//nyp=(ny+nHeight/2)%nHeight-nHeight/2;//nyp=(ny<nHeight/2) ? ny : ny-nHeight;
		nxp=(nx<=nWidth/2) ? nx : nx-nWidth;
		nyp=(ny<=nHeight/2) ? ny : ny-nHeight;
		if(0!=nxp || 0!=nyp)
		{
			tmp=(float)atan2(nyp,nxp); tsin=(float)sin(tmp); tcos=(float)cos(tmp);

			noff=ny*nWidth+nx;
			psinre[noff]=pfftim[noff]*tcos;
			psinim[noff]=-pfftre[noff]*tcos;
			tmp=pfftim[noff]*tsin;
			pfftim[noff]=-pfftre[noff]*tsin;
			pfftre[noff]=tmp;
		}
	}
	psinre[0]=psinim[0]=pfftim[0]=pfftre[0]=0;

	FFT2D(nWidth,nHeight,psinre,psinim,pcosre,pcosim,false);
	FFT2D(nWidth,nHeight,pfftre,pfftim,psinre,psinim,false);

	//Riesz Transform results calculated in the following
	for(ny=0;ny<nHeight;ny++) for(nx=0;nx<nWidth;nx++)
	{
		noff=ny*nWidth+nx;
		//Hilbert transform component
		ptOut[noff]=(float)sqrt(pow(pcosre[noff],2)+pow(psinre[noff],2));
		//Simple phase unwrapping
		nxp=nx-nWidth/2; nyp=ny-nHeight/2;
		if(0!=nyp || 0!=nxp)
		{
			tmp=(float)atan2(nyp,nxp);
			tcos=(float)cos(tmp)*pcosre[noff];
			tsin=(float)sin(tmp)*psinre[noff];
			if(abs(nyp)>abs(nxp)) {if(0>tsin) ptOut[noff]*=-1;}
			else {if(0>tcos) ptOut[noff]*=-1;}
		}
	}
	for(nx=0;nx<nsize;nx++)
	{
		if(modul==mode) ptOut[nx]=(float)sqrt(pow(ptOut[nx],2)+pow(ptIn[nx],2));//Modul
		else if(phase==mode) ptOut[nx]=(float)atan2(ptOut[nx],ptIn[nx]);//Local Phase
		else if(orientation==mode)
		{
			tsin=psinre[nx]; tcos=pcosre[nx];
			if(0!=ptOut[nx]) {tsin/=ptOut[nx]; tcos/=ptOut[nx];}
			ptOut[nx]=(float)atan2(tsin,tcos);//Orientation
		}
		else if(RTCOS==mode) ptOut[nx]=pcosre[nx];//cos part
		else if(RTSIN==mode) ptOut[nx]=psinre[nx];//sin part
	}

	delete[] pfftre; delete[] pfftim;
	delete[] pcosre; delete[] pcosim;
	delete[] psinre; delete[] psinim;
}

void CxHilbert2D::Hilbert2D(int nWidth,int nHeight,double* ptIn,double* ptOut,int mode)
{
	int nx,ny,nsize=nWidth*nHeight;

	//////////////////////////////////////////////////////////////////////
	//calculate directly the 2d transform
	/*double *pxre=new double[nsize],*pxim=new double[nsize],*pyim=new double[nsize];
	FFT2D(nWidth,nHeight,ptIn,NULL,pxre,pxim);
	for(ny=0;ny<nHeight;ny++) pxim[ny*nWidth]=pxre[ny*nWidth]=0;
	for(ny=0;ny<nHeight;ny++) for(nx=1;nx<nWidth;nx++)
	{
		int noff=ny*nWidth+nx;
		if(nx>nWidth/2) pxim[noff]*=-1;
		else pxre[noff]*=-1;
	}
	FFT2D(nWidth,nHeight,pxim,pxre,ptOut,pyim,false);
	delete[] pxre; delete[] pxim; delete[] pyim;*/

	//calculate the 2d transform via 1d Hilbert transform
	CxFFT fft;
	double *ptFTRe=new double[nWidth],*ptFTIm=new double[nWidth],*ptHTRe=new double[nWidth];
	double *ptx=ptIn,*pty=ptOut;
	for(ny=0;ny<nHeight;ny++)
	{
		fft.FFT(nWidth,ptx,NULL,ptFTRe,ptFTIm);
		if(1==nWidth%2)
		{
			for(nx=1;nx<=(nWidth-1)/2;nx++) {ptFTRe[nx]*=2; ptFTIm[nx]*=2;}
			for(nx=(nWidth+1)/2;nx<nWidth;nx++) ptFTRe[nx]=ptFTIm[nx]=0;
		}
		else
		{
			for(nx=1;nx<=nWidth/2;nx++) {ptFTRe[nx]*=2; ptFTIm[nx]*=2;}
			for(nx=nWidth/2+1;nx<nWidth;nx++) ptFTRe[nx]=ptFTIm[nx]=0;
		}
		fft.FFT(nWidth,ptFTRe,ptFTIm,ptHTRe,pty,false);
		ptx+=nWidth; pty+=nWidth;
	}
	delete[] ptFTRe; delete[] ptFTIm; delete[] ptHTRe;
	//////////////////////////////////////////////////////////////////////

	for(nx=0;nx<nsize;nx++)
	{
		if(modul==mode) ptOut[nx]=sqrt(ptIn[nx]*ptIn[nx]+ptOut[nx]*ptOut[nx]);
		else if(phase==mode) ptOut[nx]=atan2(ptOut[nx],ptIn[nx]);
	}
}

void CxHilbert2D::Riesz(int nWidth,int nHeight,double* ptIn,double* ptOut,int mode)
{
	int nx,ny,nxp,nyp,noff,nsize=nWidth*nHeight;
	double tmp,*pfftre=new double[nsize],*pfftim=new double[nsize];
	double tcos,*pcosre=new double[nsize],*pcosim=new double[nsize];
	double tsin,*psinre=new double[nsize],*psinim=new double[nsize];

	FFT2D(nWidth,nHeight,ptIn,NULL,pfftre,pfftim);

	for(ny=0;ny<nHeight;ny++) for(nx=0;nx<nWidth;nx++)
	{
		//nxp=(nx+nWidth/2)%nWidth-nWidth/2;//nxp=(nx<nWidth/2) ? nx : nx-nWidth;
		//nyp=(ny+nHeight/2)%nHeight-nHeight/2;//nyp=(ny<nHeight/2) ? ny : ny-nHeight;
		nxp=(nx<=nWidth/2) ? nx : nx-nWidth;
		nyp=(ny<=nHeight/2) ? ny : ny-nHeight;
		if(0!=nxp || 0!=nyp)
		{
			tmp=atan2(nyp,nxp); tsin=sin(tmp); tcos=cos(tmp);

			noff=ny*nWidth+nx;
			psinre[noff]=pfftim[noff]*tcos;
			psinim[noff]=-pfftre[noff]*tcos;
			tmp=pfftim[noff]*tsin;
			pfftim[noff]=-pfftre[noff]*tsin;
			pfftre[noff]=tmp;
		}
	}
	psinre[0]=psinim[0]=pfftim[0]=pfftre[0]=0;

	FFT2D(nWidth,nHeight,psinre,psinim,pcosre,pcosim,false);
	FFT2D(nWidth,nHeight,pfftre,pfftim,psinre,psinim,false);

	//Riesz Transform results calculated in the following
	for(ny=0;ny<nHeight;ny++) for(nx=0;nx<nWidth;nx++)
	{
		noff=ny*nWidth+nx;
		//Hilbert transform component
		ptOut[noff]=sqrt(pow(pcosre[noff],2)+pow(psinre[noff],2));
		//Simple phase unwrapping
		nxp=nx-nWidth/2; nyp=ny-nHeight/2;
		if(0!=nyp || 0!=nxp)
		{
			tmp=atan2(nyp,nxp);
			tcos=cos(tmp)*pcosre[noff];
			tsin=sin(tmp)*psinre[noff];
			if(abs(nyp)>abs(nxp)) {if(0>tsin) ptOut[noff]*=-1;}
			else {if(0>tcos) ptOut[noff]*=-1;}
		}
	}
	for(nx=0;nx<nsize;nx++)
	{
		if(modul==mode) ptOut[nx]=sqrt(pow(ptOut[nx],2)+pow(ptIn[nx],2));//Modul
		else if(phase==mode) ptOut[nx]=atan2(ptOut[nx],ptIn[nx]);//Local Phase
		else if(orientation==mode)
		{
			tsin=psinre[nx]; tcos=pcosre[nx];
			if(0!=ptOut[nx]) {tsin/=ptOut[nx]; tcos/=ptOut[nx];}
			ptOut[nx]=atan2(tsin,tcos);//Orientation
		}
		else if(RTCOS==mode) ptOut[nx]=pcosre[nx];//cos part
		else if(RTSIN==mode) ptOut[nx]=psinre[nx];//sin part
	}

	delete[] pfftre; delete[] pfftim;
	delete[] pcosre; delete[] pcosim;
	delete[] psinre; delete[] psinim;
}
