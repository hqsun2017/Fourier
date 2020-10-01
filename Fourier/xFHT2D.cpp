#include "xFHT2D.h"

bool CxFHT2D_s::FHT2D(int nX,int nY,float* ptHT)
{
	if(!powerof2(nY) || !powerof2(nX) ) return false;

	float* ptmp=ptHT;
	for(int y=0;y<nY;y++)
	{
		FHT(nX,ptmp);
		ptmp+=nX;
	}
	ptmp=new float[nY];
	for(int x=0;x<nX;x++)
	{
		for(y=0;y<nY;y++) ptmp[y]=ptHT[y*nX+x];
		FHT(nY,ptmp);
		for(y=0;y<nY;y++) ptHT[y*nX+x]=ptmp[y];
	}
	delete[] ptmp;

	return true;
}

bool CxFHT2D_s::Correlation(int nX,int nY,float* ptx,float* pty,float* ptc)
{
	if(false==FHT2D(nX,nY,ptx)) return false;
	if(false==FHT2D(nX,nY,pty)) return false;

	int x,y,xp,yp,ny0,ny1,n00,n11,n01,n10,ntot=nX*nY;
	for(y=0;y<nY;y++)
	{
		yp=(nY-y)%nY; ny0=y*nX; ny1=yp*nX;
		for(x=0;x<nX;x++)
		{
			xp=(nX-x)%nX;
			n00=ny0+x; n11=ny1+xp;
			n01=ny0+xp; n10=ny1+x;
			ptc[n00]=((ptx[n00]+ptx[n11])*(pty[n00]+pty[n11])
				+(ptx[n01]-ptx[n10])*(pty[n01]-pty[n10])
				+(ptx[n10]+ptx[n01])*(pty[n00]-pty[n11])
				-(ptx[n00]-ptx[n11])*(pty[n01]+pty[n10]))/4;
		}
	}
	FHT2D(nX,nY,ptc);
	for(y=0;y<ntot;y++) ptc[y]/=ntot;//////////////////////////

	return true;
}

bool CxFHT2D::FHT2D(int nX,int nY,double* ptHT)
{
	if(!powerof2(nY) || !powerof2(nX) ) return false;

	double* ptmp=ptHT;
	for(int y=0;y<nY;y++)
	{
		FHT(nX,ptmp);
		ptmp+=nX;
	}
	ptmp=new double[nY];
	for(int x=0;x<nX;x++)
	{
		for(y=0;y<nY;y++) ptmp[y]=ptHT[y*nX+x];
		FHT(nY,ptmp);
		for(y=0;y<nY;y++) ptHT[y*nX+x]=ptmp[y];
	}
	delete[] ptmp;

	return true;
}

bool CxFHT2D::Correlation(int nX,int nY,double* ptx,double* pty,double* ptc)
{
	if(false==FHT2D(nX,nY,ptx)) return false;
	if(false==FHT2D(nX,nY,pty)) return false;

	int x,y,xp,yp,ny0,ny1,n00,n11,n01,n10,ntot=nX*nY;
	for(y=0;y<nY;y++)
	{
		yp=(nY-y)%nY; ny0=y*nX; ny1=yp*nX;
		for(x=0;x<nX;x++)
		{
			xp=(nX-x)%nX;
			n00=ny0+x; n11=ny1+xp;
			n01=ny0+xp; n10=ny1+x;
			ptc[n00]=((ptx[n00]+ptx[n11])*(pty[n00]+pty[n11])
				+(ptx[n01]-ptx[n10])*(pty[n01]-pty[n10])
				+(ptx[n10]+ptx[n01])*(pty[n00]-pty[n11])
				-(ptx[n00]-ptx[n11])*(pty[n01]+pty[n10]))/4;
		}
	}
	FHT2D(nX,nY,ptc);
	for(y=0;y<ntot;y++) ptc[y]/=ntot;//////////////////////////

	return true;
}
