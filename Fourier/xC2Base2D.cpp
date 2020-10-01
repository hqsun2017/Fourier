#include "xC2Base2D.h"

void CxC2Base2D_s::LinearProcess(int nw1,int nh1,float* ptx,int nw2,int nh2,float* pty,float* ptc)
{//线性运算
	int nx,ny,nw=nw1+nw2-1,nh=nh1+nh2-1;
	float *ptxnew=new float[nw*nh],*ptynew=new float[nw*nh];

	for(nx=0;nx<nw*nh;nx++) ptxnew[nx]=ptynew[nx]=0;
	for(ny=0;ny<nh1;ny++) for(nx=0;nx<nw1;nx++)
		ptxnew[ny*nw+nx]=ptx[ny*nw1+nx];
	for(ny=0;ny<nh2;ny++) for(nx=0;nx<nw2;nx++)
		ptynew[ny*nw+nx]=pty[ny*nw2+nx];

	CircleProcess(nw,nh,ptxnew,ptynew,ptc);

	delete[] ptxnew; delete[] ptynew;
}

void CxC2Base2D::LinearProcess(int nw1,int nh1,double* ptx,int nw2,int nh2,double* pty,double* ptc)
{//线性运算
	int nx,ny,nw=nw1+nw2-1,nh=nh1+nh2-1;
	double *ptxnew=new double[nw*nh],*ptynew=new double[nw*nh];

	for(nx=0;nx<nw*nh;nx++) ptxnew[nx]=ptynew[nx]=0;
	for(ny=0;ny<nh1;ny++) for(nx=0;nx<nw1;nx++)
		ptxnew[ny*nw+nx]=ptx[ny*nw1+nx];
	for(ny=0;ny<nh2;ny++) for(nx=0;nx<nw2;nx++)
		ptynew[ny*nw+nx]=pty[ny*nw2+nx];

	CircleProcess(nw,nh,ptxnew,ptynew,ptc);

	delete[] ptxnew; delete[] ptynew;
}
