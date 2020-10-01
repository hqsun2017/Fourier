#include "xRealFFT.h"

void CxRealFFT_s::RealFFT(int n,float* x1,float* x2,float* y1Re,float* y1Im,float* y2Re,float* y2Im)
{
	FFT(n,x1,x2,y1Re,y1Im);
	int i=1,nhalf=n/2;
	y2Re[0]=y1Im[0]; y2Im[0]=0; y1Im[0]=0;
	for(;i<=nhalf;i++)
	{
		y2Re[i]=(y1Im[i]+y1Im[n-i])/2;
		y2Im[i]=(y1Re[n-i]-y1Re[i])/2;
		y1Re[i]=(y1Re[i]+y1Re[n-i])/2;
		y1Im[i]=(y1Im[i]-y1Im[n-i])/2;
	}
	for(;i<n;i++)
	{
		y2Re[i]=y2Re[n-i];
		y2Im[i]=-y2Im[n-i];
		y1Re[i]=y1Re[n-i];
		y1Im[i]=-y1Im[n-i];
	}
}

void CxRealFFT::RealFFT(int n,double* x1,double* x2,double* y1Re,double* y1Im,double* y2Re,double* y2Im)
{
	FFT(n,x1,x2,y1Re,y1Im);
	int i=1,nhalf=n/2;
	y2Re[0]=y1Im[0]; y2Im[0]=0; y1Im[0]=0;
	for(;i<=nhalf;i++)
	{
		y2Re[i]=(y1Im[i]+y1Im[n-i])/2;
		y2Im[i]=(y1Re[n-i]-y1Re[i])/2;
		y1Re[i]=(y1Re[i]+y1Re[n-i])/2;
		y1Im[i]=(y1Im[i]-y1Im[n-i])/2;
	}
	for(;i<n;i++)
	{
		y2Re[i]=y2Re[n-i];
		y2Im[i]=-y2Im[n-i];
		y1Re[i]=y1Re[n-i];
		y1Im[i]=-y1Im[n-i];
	}
}
