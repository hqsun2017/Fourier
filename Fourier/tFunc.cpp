#include "tFunc.h"

void DFT_s(int N,float* pxRe,float* pxIm,float* pyRe,float* pyIm,bool bTrans)
{
	float tPI=4.0f*(float)atan(1),tRe,tIm,tmp,te=1;
	float tw,tPhase=2*tPI/N,tcos,tsin;
	int i,j;
	if(false==bTrans) te=-1;
	for(i=0;i<N;i++)
	{
		tRe=tIm=0;
		for(j=0;j<N;j++)
		{
			tw=tPhase*i*j;
			tcos=(float)cos(tw);
			tsin=-(float)sin(tw);
			tmp=(NULL==pxIm) ? 0 : pxIm[j];
			tRe+=pxRe[j]*tcos-tmp*te*tsin;
			tIm+=tmp*te*tcos+pxRe[j]*tsin;
		}
		pyRe[i]=tRe;
		pyIm[i]=tIm;
		if(false==bTrans) {pyRe[i]/=N; pyIm[i]/=-N;}
	}
}

void DFT2D_s(int nX,int nY,float* pxRe,float* pxIm,float* pyRe,float* pyIm,bool bTrans)
{
	int nx,ny;
	float *pxre=pxRe,*pxim=pxIm,*pyre=pyRe,*pyim=pyIm;
	for(ny=0;ny<nY;ny++)
	{
		DFT_s(nX,pxre,pxim,pyre,pyim,bTrans);
		pxre+=nX;
		if(NULL!=pxIm) pxim+=nX;
		pyre+=nX;
		pyim+=nX;
	}
	pxre=new float[nY];
	pxim=new float[nY];
	pyre=new float[nY];
	pyim=new float[nY];
	for(nx=0;nx<nX;nx++)
	{
		for(ny=0;ny<nY;ny++)
		{
			pxre[ny]=pyRe[ny*nX+nx];
			pxim[ny]=pyIm[ny*nX+nx];
		}
		DFT_s(nY,pxre,pxim,pyre,pyim,bTrans);
		for(ny=0;ny<nY;ny++)
		{
			pyRe[ny*nX+nx]=pyre[ny];
			pyIm[ny*nX+nx]=pyim[ny];
		}
	}
	delete[] pxre;
	delete[] pxim;
	delete[] pyre;
	delete[] pyim;
}

void HRFT_s(int nLen,float* ptIn,float tOmega,float& fRe,float& fIm)
{
	//cos(tOmega),sin(tOmega),cos(n*tOmega) and sin(n*tOmega)
	float fcos_1,fsin_1,fcosn,fsin_n;
	//Chebyshev polynomials of Type I and II
	float fT0,fT1,fT2,fU0,fU1,fU2;

	fcos_1=(float)cos(tOmega);
	fsin_1=(float)sin(tOmega);
	fT0=1;
	fT1=fcos_1;
	fU0=1;
	fU1=2*fcos_1;

	//this takes care of n=0,1
	fRe=ptIn[0]+ptIn[1]*fT1;
	fIm=ptIn[1]*fsin_1*fU0;

	fcos_1*=2;
	for(int n=2;n<nLen;++n)
	{
		fT2=fcos_1*fT1-fT0;
		fT0=fT1;
		fT1=fT2;
		fU2=fcos_1*fU1-fU0;
		fU0=fU1;
		fU1=fU2;
		fcosn=fT2;
		fsin_n=fsin_1*fU0;
		fRe+=ptIn[n]*fcosn;
		fIm+=ptIn[n]*fsin_n;
	}
}

void Rectangular_s(int N,float* ptwin)
{
	for(int i=0;i<N;i++) ptwin[i]=1;
}

void Triangular_s(int N,float* ptwin)
{
	for(int i=0;i<N;i++) ptwin[i]=1-(float)abs(2*i+1-N)/(N-1);
}

void Hanning_s(int N,float* ptwin)
{
	float t2PI=8*(float)atan(1);
	for(int i=0;i<N;i++) ptwin[i]=0.5f-0.5f*(float)cos(t2PI*i/(N-1));
}

void Hamming_s(int N,float* ptwin)
{
	float t2PI=8*(float)atan(1);
	for(int i=0;i<N;i++) ptwin[i]=0.53836f-0.46164f*(float)cos(t2PI*i/(N-1));
}

void Gauss_s(int N,float tpara,float* ptwin)
{
	for(int i=0;i<N;i++) ptwin[i]=(float)exp(-pow((i-(N-1)/2)/(tpara*(N-1)/2),2)/2);
}

void Bartlett_Hanning_s(int N,float* ptwin)
{
	float t2PI=8*(float)atan(1);
	for(int i=0;i<N;i++) ptwin[i]=0.62f-0.38f*(float)cos(t2PI*i/(N-1))-0.24f*(float)abs(2*i+1-N)/(N-1);
}

void Blackman_s(int N,float* ptwin)
{
	float t2PI=8*(float)atan(1);
	for(int i=0;i<N;i++)
	{
		float tmp=t2PI*i/(N-1);
		ptwin[i]=(7938-9240*(float)cos(tmp)+1430*(float)cos(2*tmp))/18608;
	}
}

void Nuttall_s(int N,float* ptwin)
{
	float t2PI=8*(float)atan(1);
	for(int i=0;i<N;i++)
	{
		float tmp=t2PI*i/(N-1);
		ptwin[i]=0.355768f-0.487396f*(float)cos(tmp)+0.144232f*(float)cos(2*tmp)-0.012604f*(float)cos(3*tmp);
	}
}

void Blackman_Harris_s(int N,float* ptwin)
{
	float t2PI=8*(float)atan(1);
	for(int i=0;i<N;i++)
	{
		float tmp=t2PI*i/(N-1);
		//ptwin[i]=0.35875f-0.48829f*(float)cos(tmp)+0.14128f*(float)cos(2*tmp)-0.01168f*(float)cos(3*tmp);
		ptwin[i]=0.422323f-0.49755f*(float)cos(tmp)+0.07922f*(float)cos(2*tmp);
	}
}

void Blackman_Nuttall_s(int N,float* ptwin)
{
	float t2PI=8*(float)atan(1);
	for(int i=0;i<N;i++)
	{
		float tmp=t2PI*i/(N-1);
		ptwin[i]=0.3635819f-0.4891775f*(float)cos(tmp)+0.1363995f*(float)cos(2*tmp)-0.0106411f*(float)cos(3*tmp);
	}
}

void Welch_s(int N,float* ptwin)
{
	for(int i=0;i<N;i++) ptwin[i]=1-(float)(pow(2*i+1-N,2)/pow(N-1,2));
}

void Flat_top_s(int N,float* ptwin)
{
	float t2PI=8*(float)atan(1);
	for(int i=0;i<N;i++)
	{
		float tmp=(t2PI*i)/(N-1);
		//ptwin[i]=1-1.93*cos(tmp)+1.29*cos(2*tmp)-0.388*cos(3*tmp)+0.032*cos(4*tmp);
		ptwin[i]=0.215578948f-0.41663158f*(float)cos(tmp)+0.277263158f*(float)cos(2*tmp)-0.083578947f*(float)cos(3*tmp)+0.006947368f*(float)cos(4*tmp);
	}
}

void Bohman_s(int N,float* ptwin)
{
	float tPI=4*(float)atan(1);
	for(int i=0;i<N;i++)
	{
		float tmp=(float)abs(2*i+1-N)/(N-1);
		ptwin[i]=(1-tmp)*(float)cos(tmp*tPI)+(float)sin(tmp*tPI)/tPI;
	}
}

float BesselI0_s(float x)
{
	double denominator,numerator,z;
	if(x==0) return 1;
	else
	{
		z=x*x;
		numerator=(z*(z*(z*(z*(z*(z*(z*(z*(z*(z*(z*(z*(z*(z*0.210580722890567e-22+0.380715242345326e-19)+
			0.479440257548300e-16)+0.435125971262668e-13)+0.300931127112960e-10)+0.160224679395361e-7)+
			0.654858370096785e-5)+0.202591084143397e-2)+0.463076284721000e0)+0.754337328948189e2)+
			0.830792541809429e4)+0.571661130563785e6)+0.216415572361227e8)+0.356644482244025e9)+0.144048298227235e10);
		denominator=(z*(z*(z-0.307646912682801e4)+0.347626332405882e7)-0.144048298227235e10);
	}

	return (float)(-numerator/denominator);
}

void Kaiser_s(int N,float alpha,float* ptwin)
{
	float tPI=4*(float)atan(1);
	for(int i=0;i<N;i++)
		ptwin[i]=BesselI0_s(tPI*alpha*(float)sqrt(1-pow(2.0*i/(N-1)-1,2)))/BesselI0_s(tPI*alpha);
}

void Parzen_s(int N,float* ptwin)
{
	for(int i=0;i<N;i++)
	{
		int ntmp=abs(2*i-N+1);
		float ttmp=(float)ntmp/(N-1);
		if(ntmp<=(N-1)/2) ptwin[i]=1-6*(float)pow(ttmp,2)+6*(float)pow(ttmp,3);
		else ptwin[i]=2*(float)pow(1-ttmp,3);
	}
}

void sgn_s(int nSize,float*& psgnRe,float*& psgnIm)
{
	float pi=4*(float)atan(1);
	psgnRe=new float[nSize];
	if(0==nSize%2)
	{
		psgnIm=new float[nSize];
		for(int i=0;i<nSize;i++)
		{
			if(0==i%2) {psgnRe[i]=0; psgnIm[i]=-1.0f/(float)nSize;}
			else
			{
				float phi=2.0f*pi*(float)i/(float)nSize;
				psgnRe[i]=(2.0f*(float)sin(phi)/(1.0f-(float)cos(phi)))/(float)nSize;
				psgnIm[i]=1.0f/(float)nSize;
			}
		}
	}
	else
	{
		psgnIm=NULL; psgnRe[0]=0;
		for(int i=1;i<nSize;i++)
		{
			float phi=2.0f*pi*(float)i/(float)nSize;
			psgnRe[i]=(float)sin(phi);
			if(0==i%2) psgnRe[i]-=2.0f*(float)sin(phi/2.0);
			else psgnRe[i]+=2.0f*(float)sin(phi/2.0);
			psgnRe[i]/=(1.0f-(float)cos(phi))*(float)nSize;
		}
	}
}

void SlowDCT_s(int n,float* ptx,float* pty,bool bdir)
{
	int i,k;
	float pi=4*(float)atan(1);

	for(k=0;k<n;k++)
	{
		pty[k]=(0==k%2) ? 0.5f*ptx[n-1] : -0.5f*ptx[n-1];
		pty[k]+=0.5f*ptx[0];
		for(i=1;i<n-1;i++)
			pty[k]+=ptx[i]*(float)cos(pi*double(k*i)/double(n-1));
	}
	if(false==bdir)
		for(k=0;k<n;k++) pty[k]*=2.0f/float(n-1);
}

void SlowDST_s(int n,float* ptx,float* pty,bool bdir)
{
	int i,k;
	float pi=4*(float)atan(1);

	for(k=0;k<n;k++) pty[k]=0;
	for(k=0;k<n;k++) for(i=0;i<n;i++)
		pty[k]+=ptx[i]*(float)sin(pi*double((k+1)*(i+1))/double(n+1));
	if(false==bdir)
		for(k=0;k<n;k++) pty[k]*=2.0f/float(n+1);
}

void SlowDST0_s(int n,float* ptx,float* pty,bool bdir)
{
	int i,k;
	float pi=4*(float)atan(1);

	for(k=0;k<n;k++) pty[k]=0;
	for(k=1;k<n-1;k++) for(i=1;i<n-1;i++)
		pty[k]+=ptx[i]*(float)sin(pi*double(k*i)/double(n-1));
	if(false==bdir)
		for(k=0;k<n;k++) pty[k]*=2.0f/float(n-1);
}

//For double type in the following
void DFT(int N,double* pxRe,double* pxIm,double* pyRe,double* pyIm,bool bTrans)
{
	double tPI=4*atan(1),tRe,tIm,tmp,te=1;
	double tw,tPhase=2*tPI/N,tcos,tsin;
	int i,j;
	if(false==bTrans) te=-1;
	for(i=0;i<N;i++)
	{
		tRe=tIm=0;
		for(j=0;j<N;j++)
		{
			tw=tPhase*i*j;
			tcos=cos(tw);
			tsin=-sin(tw);
			tmp=(NULL==pxIm) ? 0 : pxIm[j];
			tRe+=pxRe[j]*tcos-tmp*te*tsin;
			tIm+=tmp*te*tcos+pxRe[j]*tsin;
		}
		pyRe[i]=tRe;
		pyIm[i]=tIm;
		if(false==bTrans) {pyRe[i]/=N; pyIm[i]/=-N;}
	}
}

void DFT2D(int nX,int nY,double* pxRe,double* pxIm,double* pyRe,double* pyIm,bool bTrans)
{
	int nx,ny;
	double *pxre=pxRe,*pxim=pxIm,*pyre=pyRe,*pyim=pyIm;
	for(ny=0;ny<nY;ny++)
	{
		DFT(nX,pxre,pxim,pyre,pyim,bTrans);
		pxre+=nX;
		if(NULL!=pxIm) pxim+=nX;
		pyre+=nX;
		pyim+=nX;
	}
	pxre=new double[nY];
	pxim=new double[nY];
	pyre=new double[nY];
	pyim=new double[nY];
	for(nx=0;nx<nX;nx++)
	{
		for(ny=0;ny<nY;ny++)
		{
			pxre[ny]=pyRe[ny*nX+nx];
			pxim[ny]=pyIm[ny*nX+nx];
		}
		DFT(nY,pxre,pxim,pyre,pyim,bTrans);
		for(ny=0;ny<nY;ny++)
		{
			pyRe[ny*nX+nx]=pyre[ny];
			pyIm[ny*nX+nx]=pyim[ny];
		}
	}
	delete[] pxre;
	delete[] pxim;
	delete[] pyre;
	delete[] pyim;
}

void HRFT(int nLen,double* ptIn,double tOmega,double& fRe,double& fIm)
{
	//cos(tOmega),sin(tOmega),cos(n*tOmega) and sin(n*tOmega)
	double fcos_1,fsin_1,fcosn,fsin_n;
	//Chebyshev polynomials of Type I and II
	double fT0,fT1,fT2,fU0,fU1,fU2;

	fcos_1=cos(tOmega);
	fsin_1=sin(tOmega);
	fT0=1;
	fT1=fcos_1;
	fU0=1;
	fU1=2*fcos_1;

	//this takes care of n=0,1
	fRe=ptIn[0]+ptIn[1]*fT1;
	fIm=ptIn[1]*fsin_1*fU0;

	fcos_1*=2;
	for(int n=2;n<nLen;++n)
	{
		fT2=fcos_1*fT1-fT0;
		fT0=fT1;
		fT1=fT2;
		fU2=fcos_1*fU1-fU0;
		fU0=fU1;
		fU1=fU2;
		fcosn=fT2;
		fsin_n=fsin_1*fU0;
		fRe+=ptIn[n]*fcosn;
		fIm+=ptIn[n]*fsin_n;
	}
}

void Rectangular(int N,double* ptwin)
{
	for(int i=0;i<N;i++) ptwin[i]=1;
}

void Triangular(int N,double* ptwin)
{
	for(int i=0;i<N;i++) ptwin[i]=1-(double)abs(2*i+1-N)/(N-1);
}

void Hanning(int N,double* ptwin)
{
	double t2PI=8*atan(1);
	for(int i=0;i<N;i++) ptwin[i]=0.5-0.5*cos(t2PI*i/(N-1));
}

void Hamming(int N,double* ptwin)
{
	double t2PI=8*atan(1);
	for(int i=0;i<N;i++) ptwin[i]=0.53836-0.46164*cos(t2PI*i/(N-1));
}

void Gauss(int N,double tpara,double* ptwin)
{
	for(int i=0;i<N;i++) ptwin[i]=exp(-pow((i-(N-1)/2)/(tpara*(N-1)/2),2)/2);
}

void Bartlett_Hanning(int N,double* ptwin)
{
	double t2PI=8*atan(1);
	for(int i=0;i<N;i++) ptwin[i]=0.62-0.38*cos(t2PI*i/(N-1))-0.24*abs(2*i+1-N)/(N-1);
}

void Blackman(int N,double* ptwin)
{
	double t2PI=8*atan(1);
	for(int i=0;i<N;i++)
	{
		double tmp=t2PI*i/(N-1);
		ptwin[i]=(7938-9240*cos(tmp)+1430*cos(2*tmp))/18608;
	}
}

void Nuttall(int N,double* ptwin)
{
	double t2PI=8*atan(1);
	for(int i=0;i<N;i++)
	{
		double tmp=t2PI*i/(N-1);
		ptwin[i]=0.355768-0.487396*cos(tmp)+0.144232*cos(2*tmp)-0.012604*cos(3*tmp);
	}
}

void Blackman_Harris(int N,double* ptwin)
{
	double t2PI=8*atan(1);
	for(int i=0;i<N;i++)
	{
		double tmp=t2PI*i/(N-1);
		//ptwin[i]=0.35875-0.48829*cos(tmp)+0.14128*cos(2*tmp)-0.01168*cos(3*tmp);
		ptwin[i]=0.422323-0.49755*cos(tmp)+0.07922*cos(2*tmp);
	}
}

void Blackman_Nuttall(int N,double* ptwin)
{
	double t2PI=8*atan(1);
	for(int i=0;i<N;i++)
	{
		double tmp=t2PI*i/(N-1);
		ptwin[i]=0.3635819-0.4891775*cos(tmp)+0.1363995*cos(2*tmp)-0.0106411*cos(3*tmp);
	}
}

void Welch(int N,double* ptwin)
{
	for(int i=0;i<N;i++) ptwin[i]=1-pow(2*i+1-N,2)/pow(N-1,2);
}

void Flat_top(int N,double* ptwin)
{
	double t2PI=8*atan(1);
	for(int i=0;i<N;i++)
	{
		double tmp=(t2PI*i)/(N-1);
		//ptwin[i]=1-1.93*cos(tmp)+1.29*cos(2*tmp)-0.388*cos(3*tmp)+0.032*cos(4*tmp);
		ptwin[i]=0.215578948-0.41663158*cos(tmp)+0.277263158*cos(2*tmp)-0.083578947*cos(3*tmp)+0.006947368*cos(4*tmp);
	}
}

void Bohman(int N,double* ptwin)
{
	double tPI=4*atan(1);
	for(int i=0;i<N;i++)
	{
		double tmp=(double)abs(2*i+1-N)/(N-1);
		ptwin[i]=(1-tmp)*cos(tmp*tPI)+sin(tmp*tPI)/tPI;
	}
}

double BesselI0(double x)
{
	double denominator,numerator,z;
	if(x==0) return 1;
	else
	{
		z=x*x;
		numerator=(z*(z*(z*(z*(z*(z*(z*(z*(z*(z*(z*(z*(z*(z*0.210580722890567e-22+0.380715242345326e-19)+
			0.479440257548300e-16)+0.435125971262668e-13)+0.300931127112960e-10)+0.160224679395361e-7)+
			0.654858370096785e-5)+0.202591084143397e-2)+0.463076284721000e0)+0.754337328948189e2)+
			0.830792541809429e4)+0.571661130563785e6)+0.216415572361227e8)+0.356644482244025e9)+0.144048298227235e10);
		denominator=(z*(z*(z-0.307646912682801e4)+0.347626332405882e7)-0.144048298227235e10);
	}

	return (-numerator/denominator);
}

void Kaiser(int N,double alpha,double* ptwin)
{
	double tPI=4*atan(1);
	for(int i=0;i<N;i++)
		ptwin[i]=BesselI0(tPI*alpha*sqrt(1-pow(2.0*i/(N-1)-1,2)))/BesselI0(tPI*alpha);
}

void Parzen(int N,double* ptwin)
{
	for(int i=0;i<N;i++)
	{
		int ntmp=abs(2*i-N+1);
		double ttmp=(double)ntmp/(N-1);
		if(ntmp<=(N-1)/2) ptwin[i]=1-6*pow(ttmp,2)+6*pow(ttmp,3);
		else ptwin[i]=2*pow(1-ttmp,3);
	}
}

void sgn(int nSize,double*& psgnRe,double*& psgnIm)
{
	double pi=4.0*atan(1.0);
	psgnRe=new double[nSize];
	if(0==nSize%2)
	{
		psgnIm=new double[nSize];
		for(int i=0;i<nSize;i++)
		{
			if(0==i%2) {psgnRe[i]=0; psgnIm[i]=-1.0/(double)nSize;}
			else
			{
				double phi=2.0*pi*(double)i/(double)nSize;
				psgnRe[i]=(2.0*sin(phi)/(1.0-cos(phi)))/(double)nSize;
				psgnIm[i]=1.0/(double)nSize;
			}
		}
	}
	else
	{
		psgnIm=NULL; psgnRe[0]=0;
		for(int i=1;i<nSize;i++)
		{
			double phi=2.0*pi*(double)i/(double)nSize;
			psgnRe[i]=sin(phi);
			if(0==i%2) psgnRe[i]-=2.0*sin(phi/2.0);
			else psgnRe[i]+=2.0*sin(phi/2.0);
			psgnRe[i]/=(1.0-cos(phi))*(double)nSize;
		}
	}
}

void SlowDCT(int n,double* ptx,double* pty,bool bdir)
{
	int i,k;
	double pi=4.0*atan(1.0);

	for(k=0;k<n;k++)
	{
		pty[k]=(0==k%2) ? 0.5*ptx[n-1] : -0.5*ptx[n-1];
		pty[k]+=0.5*ptx[0];
		for(i=1;i<n-1;i++)
			pty[k]+=ptx[i]*cos(pi*double(k*i)/double(n-1));
	}
	if(false==bdir)
		for(k=0;k<n;k++) pty[k]*=2.0/double(n-1);
}

void SlowDST(int n,double* ptx,double* pty,bool bdir)
{
	int i,k;
	double pi=4.0*atan(1.0);

	for(k=0;k<n;k++) pty[k]=0;
	for(k=0;k<n;k++) for(i=0;i<n;i++)
		pty[k]+=ptx[i]*sin(pi*double((k+1)*(i+1))/double(n+1));
	if(false==bdir)
		for(k=0;k<n;k++) pty[k]*=2.0/double(n+1);
}

void SlowDST0(int n,double* ptx,double* pty,bool bdir)
{
	int i,k;
	double pi=4.0*atan(1.0);

	for(k=0;k<n;k++) pty[k]=0;
	for(k=1;k<n-1;k++) for(i=1;i<n-1;i++)
		pty[k]+=ptx[i]*sin(pi*double(k*i)/double(n-1));
	if(false==bdir)
		for(k=0;k<n;k++) pty[k]*=2.0/double(n-1);
}
