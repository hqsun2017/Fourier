#include "x_Corr1D.h"

//循环相关
void Cx_Corr1D_s::CircleCorr1D(int nSize,float* ptx,float* pty,float* ptc)
{
	float* ptmp=new float[nSize];
	for(int i=0;i<nSize;i++) ptmp[i]=ptx[(nSize-i)%nSize];

	CircleProcess(nSize,ptmp,pty,ptc);

	delete[] ptmp;
}

//线性相关
void Cx_Corr1D_s::LinearCorr1D(int nSizx,float* ptx,int nSizy,float* pty,float* ptc)
{
	int i,nSize=nSizx+nSizy-1;
	float *ptxnew=new float[nSize],*ptynew=new float[nSize];

	for(i=0;i<nSizx;i++) ptxnew[i]=ptx[i];
	for(;i<nSize;i++) ptxnew[i]=0;
	for(i=0;i<nSizy;i++) ptynew[i]=pty[i];
	for(;i<nSize;i++) ptynew[i]=0;

	CircleCorr1D(nSize,ptxnew,ptynew,ptc);

	delete[] ptxnew; delete[] ptynew;
}

//循环相关
void Cx_Corr1D::CircleCorr1D(int nSize,double* ptx,double* pty,double* ptc)
{
	double* ptmp=new double[nSize];
	for(int i=0;i<nSize;i++) ptmp[i]=ptx[(nSize-i)%nSize];

	CircleProcess(nSize,ptmp,pty,ptc);

	delete[] ptmp;
}

//线性相关
void Cx_Corr1D::LinearCorr1D(int nSizx,double* ptx,int nSizy,double* pty,double* ptc)
{
	int i,nSize=nSizx+nSizy-1;
	double *ptxnew=new double[nSize],*ptynew=new double[nSize];

	for(i=0;i<nSizx;i++) ptxnew[i]=ptx[i];
	for(;i<nSize;i++) ptxnew[i]=0;
	for(i=0;i<nSizy;i++) ptynew[i]=pty[i];
	for(;i<nSize;i++) ptynew[i]=0;

	CircleCorr1D(nSize,ptxnew,ptynew,ptc);

	delete[] ptxnew; delete[] ptynew;
}
