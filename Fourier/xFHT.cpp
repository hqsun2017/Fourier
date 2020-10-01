#include "xFHT.h"

CxFHT_s::CxFHT_s()
{
	//tSQRT2=1.414213562373095048801688724209698;
	tSQRT2=(float)sqrt(2.0);

	float tAlpha=2*(float)atan(1);//firstly tAlpha=PI/2
	for(int i=0;i<TRIG;i++)
	{
		tCosTab[i]=(float)cos(tAlpha);
		tSinTab[i]=(float)sin(tAlpha);
		tAlpha/=2;
	}
}

bool CxFHT_s::powerof2(int n)
{
	int bits=0;
	while(n) {bits+=n&1; n>>=1;}
	return (1==bits);
}

bool CxFHT_s::FHT(int nlen,float* ptHT)
{
	if(!powerof2(nlen)) return false;

	int i,j,k,k1,k4,kx,l;
	float tx0,tx1,tx2,tx3,tx4,tx5,tx6;
	float tt,tc,ts,tc1,tc2,tc3,ts1,ts2,ts3;
	float *ptf1,*ptf2,*ptf3;

	for(k=0;(1<<k)<nlen;k++); k&=1;
	for(i=1,j=0;i<nlen;i+=1)
	{
		for(l=nlen>>1;(!((j^=l)&l));l>>=1);
		if(i>j) {tx0=ptHT[i]; ptHT[i]=ptHT[j]; ptHT[j]=tx0;}
	}

	k1=1<<k; k4=k1<<2; kx=k1>>1;
	ptf1=ptHT+k1; ptf2=ptf1+k1; ptf3=ptf2+k1;

	if(0==k)
	{
		for(i=0;i<nlen;i+=k4)
		{
			tx0=ptHT[i]+ptf1[i]; tx1=ptHT[i]-ptf1[i];
			tx2=ptf2[i]+ptf3[i]; tx3=ptf2[i]-ptf3[i];
			ptf2[i]=tx0-tx2; ptHT[i]=tx0+tx2;
			ptf3[i]=tx1-tx3; ptf1[i]=tx1+tx3;
		}
	}
	else
	{
		for(i=0,j=kx;i<nlen;i+=k4,j+=k4)
		{
			tx0=ptHT[i]-ptHT[j]; tx1=ptHT[i]+ptHT[j];
			tx2=ptf1[i]-ptf1[j]; tx3=ptf1[i]+ptf1[j];
			tx4=tx1-tx3; tx1+=tx3; tx3=tx0-tx2; tx0+=tx2;
			tx5=ptf2[i]+ptf2[j]; tx2=ptf2[i]-ptf2[j];
			tx2*=tSQRT2; ptf2[j]=tx0-tx2; ptHT[j]=tx0+tx2;
			tx2=ptf3[i]+ptf3[j]; tx0=ptf3[i]-ptf3[j];
			tx0*=tSQRT2; ptf3[j]=tx3-tx0; ptf1[j]=tx3+tx0;
			tx0=tx5-tx2; tx5+=tx2;
			ptf2[i]=tx1-tx5; ptHT[i]=tx1+tx5;
			ptf3[i]=tx4-tx0; ptf1[i]=tx4+tx0;
		}
	}

	while(k4<nlen)
	{
		k+=2; k1=1<<k; k4=k1<<2; kx=k1>>1;
		ptf1=ptHT+k1; ptf2=ptf1+k1; ptf3=ptf2+k1;
		for(i=0,j=kx;i<nlen;i+=k4,j+=k4)
		{
			tx0=ptHT[i]+ptf1[i]; tx1=ptHT[i]-ptf1[i];
			tx2=ptf2[i]+ptf3[i]; tx3=ptf2[i]-ptf3[i];
			ptf2[i]=tx0-tx2; ptHT[i]=tx0+tx2;
			ptf3[i]=tx1-tx3; ptf1[i]=tx1+tx3;
			tx1=ptHT[j]-ptf1[j]; tx0=ptHT[j]+ptf1[j];
			tx3=tSQRT2*ptf3[j]; tx2=tSQRT2*ptf2[j];
			ptf2[j]=tx0-tx2; ptHT[j]=tx0+tx2;
			ptf3[j]=tx1-tx3; ptf1[j]=tx1+tx3;
		}
		tc=tCosTab[k]; ts=tSinTab[k]; tc1=1; ts1=0;

		for(l=1;l<kx;l+=1)
		{
			tt=tc1; tc1=tt*tc-ts1*ts; ts1=tt*ts+ts1*tc;
			tc2=tc1*tc1-ts1*ts1; ts2=2*(tc1*ts1);
			tc3=tc2*tc1-ts2*ts1; ts3=tc2*ts1+ts2*tc1;

			for(i=l,j=k1-l;i<nlen;i+=k4,j+=k4)
			{
				tx0=ptf1[i]*tc2+ptf1[j]*ts2;
				tx1=ptf1[i]*ts2-ptf1[j]*tc2;
				tx2=ptf2[i]*tc1+ptf2[j]*ts1;
				tx3=ptf2[i]*ts1-ptf2[j]*tc1;
				tx4=ptf3[i]*tc3+ptf3[j]*ts3;
				tx5=ptf3[i]*ts3-ptf3[j]*tc3;
				tx6=tx2-tx4; tx4+=tx2; tx2=tx3-tx5; tx5+=tx3;
				tx3=ptHT[i]-tx0; ptf3[i]=tx3+tx2; ptf1[i]=tx3-tx2;
				tx3=ptHT[i]+tx0; ptf2[i]=tx3-tx4; ptHT[i]=tx3+tx4;
				tx3=ptHT[j]-tx1; ptf3[j]=tx3-tx5; ptf1[j]=tx3+tx5;
				tx3=ptHT[j]+tx1; ptf2[j]=tx3-tx6; ptHT[j]=tx3+tx6;
			}
		}
	}

	return true;
}

bool CxFHT_s::Spectrum(int nlen,float* ptHT)
{
	if(false==FHT(nlen,ptHT)) return false;

	int i,half=nlen/2;
	for(i=0;i<=half;i++)
		ptHT[i]=(ptHT[i]*ptHT[i]+ptHT[(nlen-i)%nlen]*ptHT[(nlen-i)%nlen])/2;
	for(;i<nlen;i++) ptHT[i]=ptHT[(nlen-i)%nlen];

	return true;
}

bool CxFHT_s::Convolution(int nlen,float* ptx,float* pty,float* ptc)
{
	if(false==FHT(nlen,ptx)) return false;
	if(false==FHT(nlen,pty)) return false;
	int i,n;
	for(i=0;i<nlen;i++)
	{
		n=(nlen-i)%nlen;
		ptc[i]=((ptx[i]+ptx[n])*pty[i]+(ptx[i]-ptx[n])*pty[n])/2;
	}
	FHT(nlen,ptc);
	for(i=0;i<nlen;i++) ptc[i]/=nlen;

	return true;
}

bool CxFHT_s::Correlation(int nlen,float* ptx,float* pty,float* ptc)
{
	if(false==FHT(nlen,ptx)) return false;
	if(false==FHT(nlen,pty)) return false;
	int i,n;
	for(i=0;i<nlen;i++)
	{
		n=(nlen-i)%nlen;
		ptc[i]=((pty[i]-pty[n])*ptx[i]+(pty[i]+pty[n])*ptx[n])/2;
	}
	FHT(nlen,ptc);
	for(i=0;i<nlen;i++) ptc[i]/=nlen;

	return true;
}

CxFHT::CxFHT()
{
	//tSQRT2=1.414213562373095048801688724209698;
	tSQRT2=sqrt(2.0);

	double tAlpha=2*atan(1.0);//firstly tAlpha=PI/2
	for(int i=0;i<TRIG;i++)
	{
		tCosTab[i]=cos(tAlpha);
		tSinTab[i]=sin(tAlpha);
		tAlpha/=2;
	}
}

bool CxFHT::powerof2(int n)
{
	int bits=0;
	while(n) {bits+=n&1; n>>=1;}
	return (1==bits);
}

bool CxFHT::FHT(int nlen,double* ptHT)
{
	if(!powerof2(nlen)) return false;

	int i,j,k,k1,k4,kx,l;
	double tx0,tx1,tx2,tx3,tx4,tx5,tx6;
	double tt,tc,ts,tc1,tc2,tc3,ts1,ts2,ts3;
	double *ptf1,*ptf2,*ptf3;

	for(k=0;(1<<k)<nlen;k++); k&=1;
	for(i=1,j=0;i<nlen;i+=1)
	{
		for(l=nlen>>1;(!((j^=l)&l));l>>=1);
		if(i>j) {tx0=ptHT[i]; ptHT[i]=ptHT[j]; ptHT[j]=tx0;}
	}

	k1=1<<k; k4=k1<<2; kx=k1>>1;
	ptf1=ptHT+k1; ptf2=ptf1+k1; ptf3=ptf2+k1;

	if(0==k)
	{
		for(i=0;i<nlen;i+=k4)
		{
			tx0=ptHT[i]+ptf1[i]; tx1=ptHT[i]-ptf1[i];
			tx2=ptf2[i]+ptf3[i]; tx3=ptf2[i]-ptf3[i];
			ptf2[i]=tx0-tx2; ptHT[i]=tx0+tx2;
			ptf3[i]=tx1-tx3; ptf1[i]=tx1+tx3;
		}
	}
	else
	{
		for(i=0,j=kx;i<nlen;i+=k4,j+=k4)
		{
			tx0=ptHT[i]-ptHT[j]; tx1=ptHT[i]+ptHT[j];
			tx2=ptf1[i]-ptf1[j]; tx3=ptf1[i]+ptf1[j];
			tx4=tx1-tx3; tx1+=tx3; tx3=tx0-tx2; tx0+=tx2;
			tx5=ptf2[i]+ptf2[j]; tx2=ptf2[i]-ptf2[j];
			tx2*=tSQRT2; ptf2[j]=tx0-tx2; ptHT[j]=tx0+tx2;
			tx2=ptf3[i]+ptf3[j]; tx0=ptf3[i]-ptf3[j];
			tx0*=tSQRT2; ptf3[j]=tx3-tx0; ptf1[j]=tx3+tx0;
			tx0=tx5-tx2; tx5+=tx2;
			ptf2[i]=tx1-tx5; ptHT[i]=tx1+tx5;
			ptf3[i]=tx4-tx0; ptf1[i]=tx4+tx0;
		}
	}

	while(k4<nlen)
	{
		k+=2; k1=1<<k; k4=k1<<2; kx=k1>>1;
		ptf1=ptHT+k1; ptf2=ptf1+k1; ptf3=ptf2+k1;
		for(i=0,j=kx;i<nlen;i+=k4,j+=k4)
		{
			tx0=ptHT[i]+ptf1[i]; tx1=ptHT[i]-ptf1[i];
			tx2=ptf2[i]+ptf3[i]; tx3=ptf2[i]-ptf3[i];
			ptf2[i]=tx0-tx2; ptHT[i]=tx0+tx2;
			ptf3[i]=tx1-tx3; ptf1[i]=tx1+tx3;
			tx1=ptHT[j]-ptf1[j]; tx0=ptHT[j]+ptf1[j];
			tx3=tSQRT2*ptf3[j]; tx2=tSQRT2*ptf2[j];
			ptf2[j]=tx0-tx2; ptHT[j]=tx0+tx2;
			ptf3[j]=tx1-tx3; ptf1[j]=tx1+tx3;
		}
		tc=tCosTab[k]; ts=tSinTab[k]; tc1=1; ts1=0;

		for(l=1;l<kx;l+=1)
		{
			tt=tc1; tc1=tt*tc-ts1*ts; ts1=tt*ts+ts1*tc;
			tc2=tc1*tc1-ts1*ts1; ts2=2*(tc1*ts1);
			tc3=tc2*tc1-ts2*ts1; ts3=tc2*ts1+ts2*tc1;

			for(i=l,j=k1-l;i<nlen;i+=k4,j+=k4)
			{
				tx0=ptf1[i]*tc2+ptf1[j]*ts2;
				tx1=ptf1[i]*ts2-ptf1[j]*tc2;
				tx2=ptf2[i]*tc1+ptf2[j]*ts1;
				tx3=ptf2[i]*ts1-ptf2[j]*tc1;
				tx4=ptf3[i]*tc3+ptf3[j]*ts3;
				tx5=ptf3[i]*ts3-ptf3[j]*tc3;
				tx6=tx2-tx4; tx4+=tx2; tx2=tx3-tx5; tx5+=tx3;
				tx3=ptHT[i]-tx0; ptf3[i]=tx3+tx2; ptf1[i]=tx3-tx2;
				tx3=ptHT[i]+tx0; ptf2[i]=tx3-tx4; ptHT[i]=tx3+tx4;
				tx3=ptHT[j]-tx1; ptf3[j]=tx3-tx5; ptf1[j]=tx3+tx5;
				tx3=ptHT[j]+tx1; ptf2[j]=tx3-tx6; ptHT[j]=tx3+tx6;
			}
		}
	}

	return true;
}

bool CxFHT::Spectrum(int nlen,double* ptHT)
{
	if(false==FHT(nlen,ptHT)) return false;

	int i,half=nlen/2;
	for(i=0;i<=half;i++)
		ptHT[i]=(ptHT[i]*ptHT[i]+ptHT[(nlen-i)%nlen]*ptHT[(nlen-i)%nlen])/2;
	for(;i<nlen;i++) ptHT[i]=ptHT[(nlen-i)%nlen];

	return true;
}

bool CxFHT::Convolution(int nlen,double* ptx,double* pty,double* ptc)
{
	if(false==FHT(nlen,ptx)) return false;
	if(false==FHT(nlen,pty)) return false;
	int i,n;
	for(i=0;i<nlen;i++)
	{
		n=(nlen-i)%nlen;
		ptc[i]=((ptx[i]+ptx[n])*pty[i]+(ptx[i]-ptx[n])*pty[n])/2;
	}
	FHT(nlen,ptc);
	for(i=0;i<nlen;i++) ptc[i]/=nlen;

	return true;
}

bool CxFHT::Correlation(int nlen,double* ptx,double* pty,double* ptc)
{
	if(false==FHT(nlen,ptx)) return false;
	if(false==FHT(nlen,pty)) return false;
	int i,n;
	for(i=0;i<nlen;i++)
	{
		n=(nlen-i)%nlen;
		ptc[i]=((pty[i]-pty[n])*ptx[i]+(pty[i]+pty[n])*ptx[n])/2;
	}
	FHT(nlen,ptc);
	for(i=0;i<nlen;i++) ptc[i]/=nlen;

	return true;
}
