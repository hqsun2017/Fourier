#include "xHilbert.h"

void CxHilbert_s::HilbertTransform(int nSize,float* ptx,float* pty,int nOption,float* ptw)
{
	if(0>nOption) return;
	switch(nOption)
	{
	case 0: Hilbert0(nSize,ptx,pty); break;
	case 1: Hilbert1(nSize,ptx,pty); break;
	case 2: Hilbert2(nSize,ptx,pty); break;
	case 3: Hilbert3(nSize,ptx,pty); break;
	case 4: Hilbert4(nSize,ptx,pty,ptw); break;
	default: break;
	}
}

//不作任何处理,直接利用Fourier变换对数据进行Hilbert变换
void CxHilbert_s::Hilbert0(int nSize,float* ptx,float* pty)
{
	int i;
	float *ptFTRe=new float[nSize],*ptFTIm=new float[nSize],*ptHTRe=new float[nSize];

	fft.FFT(nSize,ptx,NULL,ptFTRe,ptFTIm);
	if(1==nSize%2)
	{
		for(i=1;i<=(nSize-1)/2;i++) {ptFTRe[i]*=2; ptFTIm[i]*=2;}
		for(i=(nSize+1)/2;i<nSize;i++) ptFTRe[i]=ptFTIm[i]=0;
	}
	else
	{
		for(i=1;i<=nSize/2;i++) {ptFTRe[i]*=2; ptFTIm[i]*=2;}
		for(i=nSize/2+1;i<nSize;i++) ptFTRe[i]=ptFTIm[i]=0;
	}
	fft.FFT(nSize,ptFTRe,ptFTIm,ptHTRe,pty,false);

	delete[] ptFTRe; delete[] ptFTIm; delete[] ptHTRe;
}

//对称延拓数据序列后,利用Fourier变换对数据进行Hilbert变换
void CxHilbert_s::Hilbert1(int nSize,float* ptx,float* pty)
{
	int i,nSize2=nSize*2;
	float *ptxNew=new float[nSize2],*ptyNew=new float[nSize2];

	for(i=0;i<nSize;i++) ptxNew[i]=ptxNew[nSize2-1-i]=ptx[i];
	Hilbert0(nSize2,ptxNew,ptyNew);
	for(i=0;i<nSize;i++) pty[i]=ptyNew[i];

	delete[] ptxNew; delete[] ptyNew;
}

//反对称延拓数据序列后,利用Fourier变换对数据进行Hilbert变换
void CxHilbert_s::Hilbert2(int nSize,float* ptx,float* pty)
{
	int i,nSize2=nSize*2;
	float *ptxNew=new float[nSize2],*ptyNew=new float[nSize2];

	for(i=0;i<nSize;i++) ptxNew[i]=ptx[i];
	for(i=0;i<nSize;i++) ptxNew[nSize2-1-i]=-ptx[i];
	Hilbert0(nSize2,ptxNew,ptyNew);
	for(i=0;i<nSize;i++) pty[i]=ptyNew[i];

	delete[] ptxNew; delete[] ptyNew;
}

//等长补零延拓数据序列后,利用Fourier变换对数据进行Hilbert变换
void CxHilbert_s::Hilbert3(int nSize,float* ptx,float* pty)
{
	int i,nSize2=nSize*2;
	float *ptxNew=new float[nSize2],*ptyNew=new float[nSize2];

	for(i=0;i<nSize;i++) ptxNew[i]=ptx[i];
	for(i=nSize;i<nSize2;i++) ptxNew[i]=0;
	Hilbert0(nSize2,ptxNew,ptyNew);
	for(i=0;i<nSize;i++) pty[i]=ptyNew[i];

	delete[] ptxNew; delete[] ptyNew;
}

//加窗处理后,再利用Fourier变换对数据进行Hilbert变换
void CxHilbert_s::Hilbert4(int nSize,float* ptx,float* pty,float* ptw)
{
	float* ptmp=new float[nSize];

	for(int i=0;i<nSize;i++) ptmp[i]=ptx[i]*ptw[i];
	Hilbert0(nSize,ptmp,pty);

	delete[] ptmp;
}

void CxHilbert::HilbertTransform(int nSize,double* ptx,double* pty,int nOption,double* ptw)
{
	if(0>nOption) return;
	switch(nOption)
	{
	case 0: Hilbert0(nSize,ptx,pty); break;
	case 1: Hilbert1(nSize,ptx,pty); break;
	case 2: Hilbert2(nSize,ptx,pty); break;
	case 3: Hilbert3(nSize,ptx,pty); break;
	case 4: Hilbert4(nSize,ptx,pty,ptw); break;
	default: break;
	}
}

//不作任何处理,直接利用Fourier变换对数据进行Hilbert变换
void CxHilbert::Hilbert0(int nSize,double* ptx,double* pty)
{
	int i;
	double *ptFTRe=new double[nSize],*ptFTIm=new double[nSize],*ptHTRe=new double[nSize];

	fft.FFT(nSize,ptx,NULL,ptFTRe,ptFTIm);
	if(1==nSize%2)
	{
		for(i=1;i<=(nSize-1)/2;i++) {ptFTRe[i]*=2; ptFTIm[i]*=2;}
		for(i=(nSize+1)/2;i<nSize;i++) ptFTRe[i]=ptFTIm[i]=0;
	}
	else
	{
		for(i=1;i<=nSize/2;i++) {ptFTRe[i]*=2; ptFTIm[i]*=2;}
		for(i=nSize/2+1;i<nSize;i++) ptFTRe[i]=ptFTIm[i]=0;
	}
	fft.FFT(nSize,ptFTRe,ptFTIm,ptHTRe,pty,false);

	delete[] ptFTRe; delete[] ptFTIm; delete[] ptHTRe;
}

//对称延拓数据序列后,利用Fourier变换对数据进行Hilbert变换
void CxHilbert::Hilbert1(int nSize,double* ptx,double* pty)
{
	int i,nSize2=nSize*2;
	double *ptxNew=new double[nSize2],*ptyNew=new double[nSize2];

	for(i=0;i<nSize;i++) ptxNew[i]=ptxNew[nSize2-1-i]=ptx[i];
	Hilbert0(nSize2,ptxNew,ptyNew);
	for(i=0;i<nSize;i++) pty[i]=ptyNew[i];

	delete[] ptxNew; delete[] ptyNew;
}

//反对称延拓数据序列后,利用Fourier变换对数据进行Hilbert变换
void CxHilbert::Hilbert2(int nSize,double* ptx,double* pty)
{
	int i,nSize2=nSize*2;
	double *ptxNew=new double[nSize2],*ptyNew=new double[nSize2];

	for(i=0;i<nSize;i++) ptxNew[i]=ptx[i];
	for(i=0;i<nSize;i++) ptxNew[nSize2-1-i]=-ptx[i];
	Hilbert0(nSize2,ptxNew,ptyNew);
	for(i=0;i<nSize;i++) pty[i]=ptyNew[i];

	delete[] ptxNew; delete[] ptyNew;
}

//等长补零延拓数据序列后,利用Fourier变换对数据进行Hilbert变换
void CxHilbert::Hilbert3(int nSize,double* ptx,double* pty)
{
	int i,nSize2=nSize*2;
	double *ptxNew=new double[nSize2],*ptyNew=new double[nSize2];

	for(i=0;i<nSize;i++) ptxNew[i]=ptx[i];
	for(i=nSize;i<nSize2;i++) ptxNew[i]=0;
	Hilbert0(nSize2,ptxNew,ptyNew);
	for(i=0;i<nSize;i++) pty[i]=ptyNew[i];

	delete[] ptxNew; delete[] ptyNew;
}

//加窗处理后,再利用Fourier变换对数据进行Hilbert变换
void CxHilbert::Hilbert4(int nSize,double* ptx,double* pty,double* ptw)
{
	double* ptmp=new double[nSize];

	for(int i=0;i<nSize;i++) ptmp[i]=ptx[i]*ptw[i];
	Hilbert0(nSize,ptmp,pty);

	delete[] ptmp;
}
