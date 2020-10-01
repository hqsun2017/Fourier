#include "xC2Base1D.h"

void CxC2Base1D_s::LinearProcess(int nSizx,float* ptx,int nSizy,float* pty,float* ptc)
{//线性运算
	int i,nSize=nSizx+nSizy-1;
	float *ptxnew=new float[nSize],*ptynew=new float[nSize];

	for(i=0;i<nSizx;i++) ptxnew[i]=ptx[i];
	for(;i<nSize;i++) ptxnew[i]=0;
	for(i=0;i<nSizy;i++) ptynew[i]=pty[i];
	for(;i<nSize;i++) ptynew[i]=0;

	CircleProcess(nSize,ptxnew,ptynew,ptc);

	delete[] ptxnew; delete[] ptynew;
}

void CxC2Base1D::LinearProcess(int nSizx,double* ptx,int nSizy,double* pty,double* ptc)
{//线性运算
	int i,nSize=nSizx+nSizy-1;
	double *ptxnew=new double[nSize],*ptynew=new double[nSize];

	for(i=0;i<nSizx;i++) ptxnew[i]=ptx[i];
	for(;i<nSize;i++) ptxnew[i]=0;
	for(i=0;i<nSizy;i++) ptynew[i]=pty[i];
	for(;i<nSize;i++) ptynew[i]=0;

	CircleProcess(nSize,ptxnew,ptynew,ptc);

	delete[] ptxnew; delete[] ptynew;
}
