#include "xConv1D.h"

//循环卷积
void CxConv1D_s::CircleProcess(int nSize,float* ptx,float* pty,float* ptc)
{
	float re1,re2,im1,im2,*ptre=new float[nSize],*ptim=new float[nSize],*ptmp=new float[nSize];

	fft.FFT(nSize,ptx,pty,ptre,ptim);
	int i=1,nhalf=nSize/2;
	ptre[0]=ptre[0]*ptim[0]; ptim[0]=0;
	for(;i<=nhalf;i++)
	{
		re1=(ptre[i]+ptre[nSize-i])/2;
		im1=(ptim[i]-ptim[nSize-i])/2;
		re2=(ptim[i]+ptim[nSize-i])/2;
		im2=(ptre[nSize-i]-ptre[i])/2;
		ptre[i]=re1*re2-im1*im2;
		ptim[i]=re1*im2+re2*im1;
	}
	for(;i<nSize;i++)
	{
		ptre[i]=ptre[nSize-i];
		ptim[i]=-ptim[nSize-i];
	}
	fft.FFT(nSize,ptre,ptim,ptc,ptmp,false);

	delete[] ptre; delete[] ptim; delete[] ptmp;
}

//重叠相加法计算线性卷积
void CxConv1D_s::OverlapConv1D(int nSizx,float* ptx,int nSizy,float* pty,float* ptc)
{
	float *ptlong,*ptshort,*ptcind=ptc,*ptctmp;
	int i,j,nsize,ntime,nlong,nshort;

	if(nSizx>nSizy) {nlong=nSizx; ptlong=ptx; nshort=nSizy; ptshort=pty;}
	else {nlong=nSizy; ptlong=pty; nshort=nSizx; ptshort=ptx;}

	ntime=nlong/nshort; nsize=2*nshort-1;
	ptctmp=new float[nsize];
	for(i=0;i<nSizx+nSizy-1;i++) ptc[i]=0;
	for(i=0;i<ntime;i++)
	{
		LinearProcess(nshort,ptshort,nshort,ptlong,ptctmp);
		for(j=0;j<nsize;j++) ptcind[j]+=ptctmp[j];
		ptlong+=nshort; ptcind+=nshort;
	}
	nlong=nlong-ntime*nshort;
	if(0!=nlong)
	{
		LinearProcess(nshort,ptshort,nlong,ptlong,ptctmp);
		for(j=0;j<nshort+nlong-1;j++) ptcind[j]+=ptctmp[j];
	}

	delete[] ptctmp;
}

//循环卷积
void CxConv1D::CircleProcess(int nSize,double* ptx,double* pty,double* ptc)
{
	double re1,re2,im1,im2,*ptre=new double[nSize],*ptim=new double[nSize],*ptmp=new double[nSize];

	fft.FFT(nSize,ptx,pty,ptre,ptim);
	int i=1,nhalf=nSize/2;
	ptre[0]=ptre[0]*ptim[0]; ptim[0]=0;
	for(;i<=nhalf;i++)
	{
		re1=(ptre[i]+ptre[nSize-i])/2;
		im1=(ptim[i]-ptim[nSize-i])/2;
		re2=(ptim[i]+ptim[nSize-i])/2;
		im2=(ptre[nSize-i]-ptre[i])/2;
		ptre[i]=re1*re2-im1*im2;
		ptim[i]=re1*im2+re2*im1;
	}
	for(;i<nSize;i++)
	{
		ptre[i]=ptre[nSize-i];
		ptim[i]=-ptim[nSize-i];
	}
	fft.FFT(nSize,ptre,ptim,ptc,ptmp,false);

	delete[] ptre; delete[] ptim; delete[] ptmp;
}

//重叠相加法计算线性卷积
void CxConv1D::OverlapConv1D(int nSizx,double* ptx,int nSizy,double* pty,double* ptc)
{
	double *ptlong,*ptshort,*ptcind=ptc,*ptctmp;
	int i,j,nsize,ntime,nlong,nshort;

	if(nSizx>nSizy) {nlong=nSizx; ptlong=ptx; nshort=nSizy; ptshort=pty;}
	else {nlong=nSizy; ptlong=pty; nshort=nSizx; ptshort=ptx;}

	ntime=nlong/nshort; nsize=2*nshort-1;
	ptctmp=new double[nsize];
	for(i=0;i<nSizx+nSizy-1;i++) ptc[i]=0;
	for(i=0;i<ntime;i++)
	{
		LinearProcess(nshort,ptshort,nshort,ptlong,ptctmp);
		for(j=0;j<nsize;j++) ptcind[j]+=ptctmp[j];
		ptlong+=nshort; ptcind+=nshort;
	}
	nlong=nlong-ntime*nshort;
	if(0!=nlong)
	{
		LinearProcess(nshort,ptshort,nlong,ptlong,ptctmp);
		for(j=0;j<nshort+nlong-1;j++) ptcind[j]+=ptctmp[j];
	}

	delete[] ptctmp;
}
