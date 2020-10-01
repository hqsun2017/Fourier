#include "xDCT_DST.h"

void CxDCT_DST_s::EasyDCT(int n,float* ptx,float* pty,bool bdir)
{
	int i,n1=n-1;
	int n2=2*n1;
	float tc,*ptre=new float[n2],*ptim=new float[n2],*ptmp=new float[n2];

	ptmp[0]=ptx[0];
	ptmp[n1]=ptx[n1];
	for(i=1;i<n1;i++)
	{
		ptmp[i]=ptx[i];
		ptmp[n2-i]=ptx[i];
	}

	FFT(n2,ptmp,NULL,ptre,ptim);

	tc=bdir ? 2 : (float)n1;
	for(i=0;i<n;i++) pty[i]=ptre[i]/tc;

	delete[] ptre; delete[] ptim; delete[] ptmp;
}

void CxDCT_DST_s::EasyDST(int n,float* ptx,float* pty,bool bdir)
{
	int i,n1=n+1;
	int n2=2*n1;
	float tc,*ptre=new float[n2],*ptim=new float[n2],*ptmp=new float[n2];

	ptmp[0]=ptmp[n1]=0;
	for(i=1;i<n1;i++)
	{
		ptmp[i]=-ptx[i-1];
		ptmp[n2-i]=ptx[i-1];
	}

	FFT(n2,ptmp,NULL,ptre,ptim);

	tc=bdir ? 2 : (float)n1;
	for(i=0;i<n;i++) pty[i]=ptim[i+1]/tc;

	delete[] ptre; delete[] ptim; delete[] ptmp;
}

void CxDCT_DST_s::EasyDST0(int n,float* ptx,float* pty,bool bdir)
{//!!!!!!the data at the begining and end of array are 0
	int i,n1=n-1;
	int n2=2*n1;
	float tc,*ptre=new float[n2],*ptim=new float[n2],*ptmp=new float[n2];

	ptmp[0]=ptmp[n1]=0;
	for(i=1;i<n1;i++)
	{
		ptmp[i]=-ptx[i];
		ptmp[n2-i]=ptx[i];
	}

	FFT(n2,ptmp,NULL,ptre,ptim);

	tc=bdir ? 2 : (float)n1;
	for(i=0;i<n;i++) pty[i]=ptim[i]/tc;

	delete[] ptre; delete[] ptim; delete[] ptmp;
}

void CxDCT_DST_s::FastDCT(int n,float* ptx1,float* ptx2,float* pty1,float* pty2,bool bdir)
{
	int i,n1=n-1;
	int n2=2*n1;
	float *ptre1=new float[n2],*ptim1=new float[n2];
	float tc,*ptre2=new float[n2],*ptim2=new float[n2];

	ptre1[0]=ptx1[0]; ptre1[n1]=ptx1[n1];
	ptim1[0]=ptx2[0]; ptim1[n1]=ptx2[n1];
	for(i=1;i<n1;i++)
	{
		ptre1[i]=ptx1[i]; ptre1[n2-i]=ptx1[i];
		ptim1[i]=ptx2[i]; ptim1[n2-i]=ptx2[i];
	}

	FFT(n2,ptre1,ptim1,ptre2,ptim2);

	tc=bdir ? 2 : (float)n1;
	for(i=0;i<n;i++)
	{
		pty1[i]=ptre2[i]/tc;
		pty2[i]=ptim2[i]/tc;
	}

	delete[] ptre1; delete[] ptim1;
	delete[] ptre2; delete[] ptim2;
}

void CxDCT_DST_s::FastDST(int n,float* ptx1,float* ptx2,float* pty1,float* pty2,bool bdir)
{
	int i,n1=n+1;
	int n2=2*n1;
	float *ptre1=new float[n2],*ptim1=new float[n2];
	float tc,*ptre2=new float[n2],*ptim2=new float[n2];

	ptre1[0]=ptre1[n1]=0;
	ptim1[0]=ptim1[n1]=0;
	for(i=1;i<n1;i++)
	{
		ptre1[i]=-ptx1[i-1]; ptre1[n2-i]=ptx1[i-1];
		ptim1[i]=-ptx2[i-1]; ptim1[n2-i]=ptx2[i-1];
	}

	FFT(n2,ptre1,ptim1,ptre2,ptim2);

	tc=bdir ? 2 : (float)n1;
	for(i=0;i<n;i++)
	{
		pty1[i]=ptim2[i+1]/tc;
		pty2[i]=-ptre2[i+1]/tc;
	}

	delete[] ptre1; delete[] ptim1;
	delete[] ptre2; delete[] ptim2;
}

void CxDCT_DST_s::FastDST0(int n,float* ptx1,float* ptx2,float* pty1,float* pty2,bool bdir)
{//!!!!!!the data at the begining and end of array are 0
	int i,n1=n-1;
	int n2=2*n1;
	float *ptre1=new float[n2],*ptim1=new float[n2];
	float tc,*ptre2=new float[n2],*ptim2=new float[n2];

	ptre1[0]=ptre1[n1]=0;
	ptim1[0]=ptim1[n1]=0;
	for(i=1;i<n1;i++)
	{
		ptre1[i]=-ptx1[i]; ptre1[n2-i]=ptx1[i];
		ptim1[i]=-ptx2[i]; ptim1[n2-i]=ptx2[i];
	}

	FFT(n2,ptre1,ptim1,ptre2,ptim2);

	tc=bdir ? 2 : (float)n1;
	for(i=0;i<n;i++)
	{
		pty1[i]=ptim2[i]/tc;
		pty2[i]=-ptre2[i]/tc;
	}

	delete[] ptre1; delete[] ptim1;
	delete[] ptre2; delete[] ptim2;
}

void CxDCT_DST::EasyDCT(int n,double* ptx,double* pty,bool bdir)
{
	int i,n1=n-1;
	int n2=2*n1;
	double tc,*ptre=new double[n2],*ptim=new double[n2],*ptmp=new double[n2];

	ptmp[0]=ptx[0];
	ptmp[n1]=ptx[n1];
	for(i=1;i<n1;i++)
	{
		ptmp[i]=ptx[i];
		ptmp[n2-i]=ptx[i];
	}

	FFT(n2,ptmp,NULL,ptre,ptim);

	tc=bdir ? 2 : n1;
	for(i=0;i<n;i++) pty[i]=ptre[i]/tc;

	delete[] ptre; delete[] ptim; delete[] ptmp;
}

void CxDCT_DST::EasyDST(int n,double* ptx,double* pty,bool bdir)
{
	int i,n1=n+1;
	int n2=2*n1;
	double tc,*ptre=new double[n2],*ptim=new double[n2],*ptmp=new double[n2];

	ptmp[0]=ptmp[n1]=0;
	for(i=1;i<n1;i++)
	{
		ptmp[i]=-ptx[i-1];
		ptmp[n2-i]=ptx[i-1];
	}

	FFT(n2,ptmp,NULL,ptre,ptim);

	tc=bdir ? 2 : n1;
	for(i=0;i<n;i++) pty[i]=ptim[i+1]/tc;

	delete[] ptre; delete[] ptim; delete[] ptmp;
}

void CxDCT_DST::EasyDST0(int n,double* ptx,double* pty,bool bdir)
{//!!!!!!the data at the begining and end of array are 0
	int i,n1=n-1;
	int n2=2*n1;
	double tc,*ptre=new double[n2],*ptim=new double[n2],*ptmp=new double[n2];

	ptmp[0]=ptmp[n1]=0;
	for(i=1;i<n1;i++)
	{
		ptmp[i]=-ptx[i];
		ptmp[n2-i]=ptx[i];
	}

	FFT(n2,ptmp,NULL,ptre,ptim);

	tc=bdir ? 2 : n1;
	for(i=0;i<n;i++) pty[i]=ptim[i]/tc;

	delete[] ptre; delete[] ptim; delete[] ptmp;
}

void CxDCT_DST::FastDCT(int n,double* ptx1,double* ptx2,double* pty1,double* pty2,bool bdir)
{
	int i,n1=n-1;
	int n2=2*n1;
	double *ptre1=new double[n2],*ptim1=new double[n2];
	double tc,*ptre2=new double[n2],*ptim2=new double[n2];

	ptre1[0]=ptx1[0]; ptre1[n1]=ptx1[n1];
	ptim1[0]=ptx2[0]; ptim1[n1]=ptx2[n1];
	for(i=1;i<n1;i++)
	{
		ptre1[i]=ptx1[i]; ptre1[n2-i]=ptx1[i];
		ptim1[i]=ptx2[i]; ptim1[n2-i]=ptx2[i];
	}

	FFT(n2,ptre1,ptim1,ptre2,ptim2);

	tc=bdir ? 2 : n1;
	for(i=0;i<n;i++)
	{
		pty1[i]=ptre2[i]/tc;
		pty2[i]=ptim2[i]/tc;
	}

	delete[] ptre1; delete[] ptim1;
	delete[] ptre2; delete[] ptim2;
}

void CxDCT_DST::FastDST(int n,double* ptx1,double* ptx2,double* pty1,double* pty2,bool bdir)
{
	int i,n1=n+1;
	int n2=2*n1;
	double *ptre1=new double[n2],*ptim1=new double[n2];
	double tc,*ptre2=new double[n2],*ptim2=new double[n2];

	ptre1[0]=ptre1[n1]=0;
	ptim1[0]=ptim1[n1]=0;
	for(i=1;i<n1;i++)
	{
		ptre1[i]=-ptx1[i-1]; ptre1[n2-i]=ptx1[i-1];
		ptim1[i]=-ptx2[i-1]; ptim1[n2-i]=ptx2[i-1];
	}

	FFT(n2,ptre1,ptim1,ptre2,ptim2);

	tc=bdir ? 2 : n1;
	for(i=0;i<n;i++)
	{
		pty1[i]=ptim2[i+1]/tc;
		pty2[i]=-ptre2[i+1]/tc;
	}

	delete[] ptre1; delete[] ptim1;
	delete[] ptre2; delete[] ptim2;
}

void CxDCT_DST::FastDST0(int n,double* ptx1,double* ptx2,double* pty1,double* pty2,bool bdir)
{//!!!!!!the data at the begining and end of array are 0
	int i,n1=n-1;
	int n2=2*n1;
	double *ptre1=new double[n2],*ptim1=new double[n2];
	double tc,*ptre2=new double[n2],*ptim2=new double[n2];

	ptre1[0]=ptre1[n1]=0;
	ptim1[0]=ptim1[n1]=0;
	for(i=1;i<n1;i++)
	{
		ptre1[i]=-ptx1[i]; ptre1[n2-i]=ptx1[i];
		ptim1[i]=-ptx2[i]; ptim1[n2-i]=ptx2[i];
	}

	FFT(n2,ptre1,ptim1,ptre2,ptim2);

	tc=bdir ? 2 : n1;
	for(i=0;i<n;i++)
	{
		pty1[i]=ptim2[i]/tc;
		pty2[i]=-ptre2[i]/tc;
	}

	delete[] ptre1; delete[] ptim1;
	delete[] ptre2; delete[] ptim2;
}
