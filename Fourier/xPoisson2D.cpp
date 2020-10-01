#include "xPoisson2D.h"

CxPoisson2D_s::CxPoisson2D_s()
{
	m_nWidth=m_nHeight=0; m_tdx=m_tdy=0;
	m_txSuperior=m_tySuperior=1; m_txInferior=0;
	m_ptSource=m_ptSolver=NULL;
	m_tAlfaL=m_tBetaL=m_tAlfaH=m_tBetaH=0;
	m_ptBondL=m_ptBondH=NULL;
}

bool CxPoisson2D_s::Check1()
{
	if(0==m_nWidth || 0==m_nHeight) return false;
	if(0==m_tdx || 0==m_tdy) return false;
	if(NULL==m_ptSolver) return false;
	if((0!=m_tAlfaL || 0!=m_tBetaL) && NULL==m_ptBondL) return false;
	if((0!=m_tAlfaH || 0!=m_tBetaH) && NULL==m_ptBondH) return false;
	return true;
}

bool CxPoisson2D_s::Check2()
{
	if(fabs(m_tySuperior-m_tdy*(m_nHeight-1))>0.5*m_tdy) return false;
	if(fabs((m_txSuperior-m_txInferior)-m_tdx*(m_nWidth-1))>0.5*m_tdx) return false;
	return true;
}

//Function to invert tridiagonal matrix equation.Matrix is nlen by nlen.
//Left,centre and right diagonal elements is stored in pta,ptb,ptc.
//Right-hand side is stored in array ptw.Solution is written to array ptu.
//pta,ptb,ptc,ptw,ptu are of extent nlen+2 with redundant 0 and nlen+1 elements.
void CxPoisson2D_s::Tridiagonal(int ndim,float* pta,float* ptb,float* ptc,float* ptw,float*& ptu)
{
	int i,nlen=ndim-2;
	float *ptx=new float[nlen],*pty=new float[nlen];

	//Scan up diagonal from i=nlen to 1
	ptx[nlen-1]=-pta[nlen]/ptb[nlen];
	pty[nlen-1]=ptw[nlen]/ptb[nlen];
	for(i=nlen-2;i>0;i--)
	{
		ptx[i]=-pta[i+1]/(ptb[i+1]+ptc[i+1]*ptx[i+1]);
		pty[i]=(ptw[i+1]-ptc[i+1]*pty[i+1])/(ptb[i+1]+ptc[i+1]*ptx[i+1]);
	}
	ptx[0]=0; pty[0]=(ptw[1]-ptc[1]*pty[1])/(ptb[1]+ptc[1]*ptx[1]);

	//Scan down diagonal from i=1 to nlen
	ptu[1]=pty[0];
	for(i=1;i<nlen;i++) ptu[i+1]=ptx[i]*ptu[i]+pty[i];

	delete[] ptx; delete[] pty;
}

bool CxPoisson2D_s::Periodic()
{
	if(false==Check1()) return false;

	int nx,ny,nxid,nyid,noff,nsize=m_nWidth*m_nHeight;
	float txx,tyy,t2pi=8*(float)atan(1);
	float *ptre=new float[nsize],*ptim=new float[nsize],*ptmp=new float[nsize];

	//Since the condition is periodic,the mean value of m_ptSource should be zero.
	fft2d.FFT2D(m_nWidth,m_nHeight,m_ptSource,NULL,ptre,ptim);
	for(ny=0;ny<m_nHeight;ny++)
	{
		nyid=ny>m_nHeight/2 ? ny-m_nHeight : ny;
		noff=ny*m_nWidth;
		for(nx=0;nx<m_nWidth;nx++)
		{
			if(0==nx && 0==ny) continue;
			nxid=nx>m_nWidth/2 ? nx-m_nWidth : nx;
			txx=2*((float)cos(t2pi*nxid/m_nWidth)-1)/(float)pow(m_tdx,2.0);
			tyy=2*((float)cos(t2pi*nyid/m_nHeight)-1)/(float)pow(m_tdy,2.0);
			ptre[noff+nx]/=txx+tyy; ptim[noff+nx]/=txx+tyy;
		}
	}
	fft2d.FFT2D(m_nWidth,m_nHeight,ptre,ptim,m_ptSolver,ptmp,false);
	delete[] ptre; delete[] ptim; delete[] ptmp;

	return true;
}

bool CxPoisson2D_s::EasyDirichlet()
{
	if(false==Check1() || false==Check2()) return false;

	int nx,ny,nwidth=m_nWidth-2,nheight=m_nHeight-1;
	float *ptV=new float[m_nWidth*m_nHeight],*ptU=new float[m_nWidth*m_nHeight];
	float *ptgamaL=new float[m_nHeight],*ptgamaH=new float[m_nHeight];
	float *ptIn=new float[m_nHeight],*ptOut=new float[m_nHeight];
	float *pta=new float[m_nWidth],*ptb=new float[m_nWidth],*ptc=new float[m_nWidth];
	float tpi=4*(float)atan(1),*ptw=new float[m_nWidth],*ptu=new float[m_nWidth];
	float tkappa=tpi*m_tdx/m_tySuperior;

	//Inverse Fourier sine transform boundary conditions
	tran.EasyDST0(m_nHeight,m_ptBondL,ptgamaL,false);
	tran.EasyDST0(m_nHeight,m_ptBondH,ptgamaH,false);

	for(nx=1;nx<=nwidth;nx++)
	{
		if(NULL!=m_ptSource)//Inverse Fourier sine transform source term
		{
			for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=m_ptSource[ny*m_nWidth+nx];
			tran.EasyDST0(m_nHeight,ptIn,ptOut,false);
			for(ny=0;ny<m_nHeight;ny++) ptV[ny*m_nWidth+nx]=ptOut[ny];
		}
		else for(ny=0;ny<m_nHeight;ny++) ptV[ny*m_nWidth+nx]=0;
	}

	for(ny=1;ny<nheight;ny++)//Solve tridiagonal matrix equations
	{
		//Initialize tridiagonal matrix
		for(nx=2;nx<=nwidth;nx++) pta[nx]=1;
		for(nx=1;nx<=nwidth;nx++) ptb[nx]=-2-(float)pow(ny*tkappa,2.0);
		ptb[1]-=m_tBetaL/(m_tAlfaL*m_tdx-m_tBetaL);
		ptb[nwidth]+=m_tBetaH/(m_tAlfaH*m_tdx+m_tBetaH);
		for(nx=1;nx<nwidth;nx++) ptc[nx]=1;

		//Initialize right-hand side vector
		for(nx=1;nx<=nwidth;nx++) ptw[nx]=ptV[ny*m_nWidth+nx]*m_tdx*m_tdx;
		ptw[1]-=ptgamaL[ny]*m_tdx/(m_tAlfaL*m_tdx-m_tBetaL);
		ptw[nwidth]-=ptgamaH[ny]*m_tdx/(m_tAlfaH*m_tdx+m_tBetaH);

		//Invert tridiagonal matrix equation
		Tridiagonal(m_nWidth,pta,ptb,ptc,ptw,ptu);
		for(nx=1;nx<=nwidth;nx++) ptU[ny*m_nWidth+nx]=ptu[nx];
	}

	for(nx=1;nx<=nwidth;nx++) ptU[nx]=ptU[nheight*m_nWidth+nx]=0;

	for(nx=1;nx<=nwidth;nx++)//Reconstruct solution via Fourier sine transform
	{
		for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=ptU[ny*m_nWidth+nx];
		tran.EasyDST0(m_nHeight,ptIn,ptOut);
		for(ny=0;ny<m_nHeight;ny++) m_ptSolver[ny*m_nWidth+nx]=ptOut[ny];
	}

	for(ny=0;ny<m_nHeight;ny++)//Calculate nx=0 and nx=m_nWidth-1 values
	{
		m_ptSolver[ny*m_nWidth]=(m_ptBondL[ny]*m_tdx-m_tBetaL*m_ptSolver[ny*m_nWidth+1])/(m_tAlfaL*m_tdx-m_tBetaL);
		m_ptSolver[(ny+1)*m_nWidth-1]=(m_ptBondH[ny]*m_tdx+m_tBetaH*m_ptSolver[ny*m_nWidth+nwidth])/(m_tAlfaH*m_tdx+m_tBetaH);
	}

	delete[] ptV; delete[] ptU;
	delete[] ptgamaL; delete[] ptgamaH;
	delete[] ptIn; delete[] ptOut;
	delete[] pta; delete[] ptb; delete[] ptc;
	delete[] ptw; delete[] ptu;

	return true;
}

bool CxPoisson2D_s::EasyNeumann()
{
	if(false==Check1() || false==Check2()) return false;

	int nx,ny,nwidth=m_nWidth-2,nheight=m_nHeight-1;
	float *ptV=new float[m_nWidth*m_nHeight],*ptU=new float[m_nWidth*m_nHeight];
	float *ptgamaL=new float[m_nHeight],*ptgamaH=new float[m_nHeight];
	float *ptIn=new float[m_nHeight],*ptOut=new float[m_nHeight];
	float *pta=new float[m_nWidth],*ptb=new float[m_nWidth],*ptc=new float[m_nWidth];
	float tpi=4*(float)atan(1),*ptw=new float[m_nWidth],*ptu=new float[m_nWidth];
	float tkappa=tpi*m_tdx/m_tySuperior;

	//Inverse Fourier cosine transform boundary conditions
	tran.EasyDCT(m_nHeight,m_ptBondL,ptgamaL,false);
	tran.EasyDCT(m_nHeight,m_ptBondH,ptgamaH,false);

	for(nx=1;nx<=nwidth;nx++)
	{
		if(NULL!=m_ptSource)//Inverse Fourier cosine transform source term
		{
			for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=m_ptSource[ny*m_nWidth+nx];
			tran.EasyDCT(m_nHeight,ptIn,ptOut,false);
			for(ny=0;ny<m_nHeight;ny++) ptV[ny*m_nWidth+nx]=ptOut[ny];
		}
		else for(ny=0;ny<m_nHeight;ny++) ptV[ny*m_nWidth+nx]=0;
	}

	for(ny=0;ny<m_nHeight;ny++)//Solve tridiagonal matrix equations
	{	
		//Initialize tridiagonal matrix
		for(nx=2;nx<=nwidth;nx++) pta[nx]=1;
		for(nx=1;nx<=nwidth;nx++) ptb[nx]=-2-(float)pow(ny*tkappa,2.0);
		ptb[1]-=m_tBetaL/(m_tAlfaL*m_tdx-m_tBetaL);
		ptb[nwidth]+=m_tBetaH/(m_tAlfaH*m_tdx+m_tBetaH);
		for(nx=1;nx<nwidth;nx++) ptc[nx]=1;
		
		//Initialize right-hand side vector
		for(nx=1;nx<=nwidth;nx++) ptw[nx]=ptV[ny*m_nWidth+nx]*m_tdx*m_tdx;
		ptw[1]-=ptgamaL[ny]*m_tdx/(m_tAlfaL*m_tdx-m_tBetaL);
		ptw[nwidth]-=ptgamaH[ny]*m_tdx/(m_tAlfaH*m_tdx+m_tBetaH);
		
		//Invert tridiagonal matrix equation
		Tridiagonal(m_nWidth,pta,ptb,ptc,ptw,ptu);
		for(nx=1;nx<=nwidth;nx++) ptU[ny*m_nWidth+nx]=ptu[nx];
	}

	for(nx=1;nx<=nwidth;nx++)//Reconstruct solution via Fourier cosine transform
	{
		for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=ptU[ny*m_nWidth+nx];
		tran.EasyDCT(m_nHeight,ptIn,ptOut);
		for(ny=0;ny<m_nHeight;ny++) m_ptSolver[ny*m_nWidth+nx]=ptOut[ny];
	}

	for(ny=0;ny<m_nHeight;ny++)//Calculate nx=0 and nx=m_nWidth-1 values
	{
		m_ptSolver[ny*m_nWidth]=(m_ptBondL[ny]*m_tdx-m_tBetaL*m_ptSolver[ny*m_nWidth+1])/(m_tAlfaL*m_tdx-m_tBetaL);
		m_ptSolver[(ny+1)*m_nWidth-1]=(m_ptBondH[ny]*m_tdx+m_tBetaH*m_ptSolver[ny*m_nWidth+nwidth])/(m_tAlfaH*m_tdx+m_tBetaH);
	}

	delete[] ptV; delete[] ptU;
	delete[] ptgamaL; delete[] ptgamaH;
	delete[] ptIn; delete[] ptOut;
	delete[] pta; delete[] ptb; delete[] ptc;
	delete[] ptw; delete[] ptu;

	return true;
}

bool CxPoisson2D_s::FastDirichlet()
{
	if(false==Check1() || false==Check2()) return false;

	int nx,ny,nwidth=m_nWidth-2,nheight=m_nHeight-1;
	float *ptV=new float[m_nWidth*m_nHeight],*ptU=new float[m_nWidth*m_nHeight];
	float *ptgamaL=new float[m_nHeight],*ptgamaH=new float[m_nHeight];
	float *ptIn=new float[m_nHeight],*ptOut=new float[m_nHeight];
	float *pta=new float[m_nWidth],*ptb=new float[m_nWidth],*ptc=new float[m_nWidth];
	float tpi=4*(float)atan(1),*ptw=new float[m_nWidth],*ptu=new float[m_nWidth];
	float tkappa=tpi*m_tdx/m_tySuperior;

	int noff,nh2=2*nheight;
	float *ptre0=new float[nh2],*ptim0=new float[nh2];
	float *ptre1=new float[nh2],*ptim1=new float[nh2];

	//Inverse Fourier sine transform boundary conditions
	ptre0[0]=ptre0[nheight]=0;
	ptim0[0]=ptim0[nheight]=0;
	for(ny=1;ny<nheight;ny++)
	{
		ptre0[ny]=-m_ptBondL[ny];
		ptre0[nh2-ny]=m_ptBondL[ny];
		ptim0[ny]=-m_ptBondH[ny];
		ptim0[nh2-ny]=m_ptBondH[ny];
	}
	tran.FFT(nh2,ptre0,ptim0,ptre1,ptim1);
	for(ny=0;ny<m_nHeight;ny++)
	{
		ptgamaL[ny]=ptim1[ny]/nheight;
		ptgamaH[ny]=-ptre1[ny]/nheight;
	}

	if(NULL!=m_ptSource)//Inverse Fourier sine transform source term
	{
		for(nx=1;nx<nwidth;nx+=2)
		{
			ptre0[0]=ptre0[nheight]=0;
			ptim0[0]=ptim0[nheight]=0;
			for(ny=1;ny<nheight;ny++)
			{
				noff=ny*m_nWidth+nx;
				ptre0[ny]=-m_ptSource[noff];
				ptre0[nh2-ny]=m_ptSource[noff];
				ptim0[ny]=-m_ptSource[noff+1];
				ptim0[nh2-ny]=m_ptSource[noff+1];
			}
			tran.FFT(nh2,ptre0,ptim0,ptre1,ptim1);
			for(ny=0;ny<m_nHeight;ny++)
			{
				noff=ny*m_nWidth+nx;
				ptV[noff]=ptim1[ny]/nheight;
				ptV[noff+1]=-ptre1[ny]/nheight;
			}
		}
		if(0!=nwidth%2)
		{
			for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=m_ptSource[ny*m_nWidth+nwidth];
			tran.EasyDST0(m_nHeight,ptIn,ptOut,false);
			for(ny=0;ny<m_nHeight;ny++) ptV[ny*m_nWidth+nwidth]=ptOut[ny];
		}
	}
	else for(nx=1;nx<=nwidth;nx++) for(ny=0;ny<m_nHeight;ny++) ptV[ny*m_nWidth+nx]=0;

	//Solve tridiagonal matrix equations
	for(ny=1;ny<nheight;ny++)
	{
		//Initialize tridiagonal matrix
		for(nx=2;nx<=nwidth;nx++) pta[nx]=1;
		for(nx=1;nx<=nwidth;nx++) ptb[nx]=-2-(float)pow(ny*tkappa,2.0);
		ptb[1]-=m_tBetaL/(m_tAlfaL*m_tdx-m_tBetaL);
		ptb[nwidth]+=m_tBetaH/(m_tAlfaH*m_tdx+m_tBetaH);
		for(nx=1;nx<nwidth;nx++) ptc[nx]=1;

		//Initialize right-hand side vector
		for(nx=1;nx<=nwidth;nx++) ptw[nx]=ptV[ny*m_nWidth+nx]*m_tdx*m_tdx;
		ptw[1]-=ptgamaL[ny]*m_tdx/(m_tAlfaL*m_tdx-m_tBetaL);
		ptw[nwidth]-=ptgamaH[ny]*m_tdx/(m_tAlfaH*m_tdx+m_tBetaH);

		//Invert tridiagonal matrix equation
		Tridiagonal(m_nWidth,pta,ptb,ptc,ptw,ptu);
		for(nx=1;nx<=nwidth;nx++) ptU[ny*m_nWidth+nx]=ptu[nx];
	}
	for(nx=1;nx<=nwidth;nx++) ptU[nx]=ptU[nheight*m_nWidth+nx]=0;

	//Reconstruct solution via Fourier sine transform
	for(nx=1;nx<nwidth;nx+=2)
	{
		ptre0[0]=ptre0[nheight]=0;
		ptim0[0]=ptim0[nheight]=0;
		for(ny=1;ny<nheight;ny++)
		{
			noff=ny*m_nWidth+nx;
			ptre0[ny]=-ptU[noff];
			ptre0[nh2-ny]=ptU[noff];
			ptim0[ny]=-ptU[noff+1];
			ptim0[nh2-ny]=ptU[noff+1];
		}
		tran.FFT(nh2,ptre0,ptim0,ptre1,ptim1);
		for(ny=0;ny<m_nHeight;ny++)
		{
			noff=ny*m_nWidth+nx;
			m_ptSolver[noff]=ptim1[ny]/2;
			m_ptSolver[noff+1]=-ptre1[ny]/2;
		}
	}
	if(0!=nwidth%2)
	{
		for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=ptU[ny*m_nWidth+nwidth];
		tran.EasyDST0(m_nHeight,ptIn,ptOut);
		for(ny=0;ny<m_nHeight;ny++) m_ptSolver[ny*m_nWidth+nwidth]=ptOut[ny];
	}

	//Calculate nx=0 and nx=m_nWidth-1 values
	for(ny=0;ny<m_nHeight;ny++)
	{
		m_ptSolver[ny*m_nWidth]=(m_ptBondL[ny]*m_tdx-m_tBetaL*m_ptSolver[ny*m_nWidth+1])/(m_tAlfaL*m_tdx-m_tBetaL);
		m_ptSolver[(ny+1)*m_nWidth-1]=(m_ptBondH[ny]*m_tdx+m_tBetaH*m_ptSolver[ny*m_nWidth+nwidth])/(m_tAlfaH*m_tdx+m_tBetaH);
	}

 	delete[] ptre0; delete[] ptim0;
	delete[] ptre1; delete[] ptim1;
	delete[] ptV; delete[] ptU;
	delete[] ptgamaL; delete[] ptgamaH;
	delete[] ptIn; delete[] ptOut;
	delete[] pta; delete[] ptb; delete[] ptc;
	delete[] ptw; delete[] ptu;

	return true;
}

bool CxPoisson2D_s::FastNeumann()
{
	if(false==Check1() || false==Check2()) return false;

	int nx,ny,nwidth=m_nWidth-2,nheight=m_nHeight-1;
	float *ptV=new float[m_nWidth*m_nHeight],*ptU=new float[m_nWidth*m_nHeight];
	float *ptgamaL=new float[m_nHeight],*ptgamaH=new float[m_nHeight];
	float *ptIn=new float[m_nHeight],*ptOut=new float[m_nHeight];
	float *pta=new float[m_nWidth],*ptb=new float[m_nWidth],*ptc=new float[m_nWidth];
	float tpi=4*(float)atan(1),*ptw=new float[m_nWidth],*ptu=new float[m_nWidth];
	float tkappa=tpi*m_tdx/m_tySuperior;

	//Inverse Fourier cosine transform boundary conditions
	int noff,nh2=2*nheight;
	float *ptre0=new float[nh2],*ptim0=new float[nh2];
	float *ptre1=new float[nh2],*ptim1=new float[nh2];

	for(ny=0;ny<nheight;ny++)
	{
		ptre0[ny]=m_ptBondL[ny];
		ptre0[nh2-1-ny]=m_ptBondL[ny+1];
		ptim0[ny]=m_ptBondH[ny];
		ptim0[nh2-1-ny]=m_ptBondH[ny+1];
	}
	tran.FFT(nh2,ptre0,ptim0,ptre1,ptim1);
	for(ny=0;ny<m_nHeight;ny++)
	{
		ptgamaL[ny]=ptre1[ny]/nheight;
		ptgamaH[ny]=ptim1[ny]/nheight;
	}

	if(NULL!=m_ptSource)//Inverse Fourier cosine transform source term
	{
		for(nx=1;nx<nwidth;nx+=2)
		{
			ptre0[0]=m_ptSource[nx];
			ptim0[0]=m_ptSource[nx+1];
			noff=nheight*m_nWidth+nx;
			ptre0[nheight]=m_ptSource[noff];
			ptim0[nheight]=m_ptSource[noff+1];
			for(ny=1;ny<nheight;ny++)
			{
				noff=ny*m_nWidth+nx;
				ptre0[ny]=m_ptSource[noff];
				ptre0[nh2-ny]=m_ptSource[noff];
				ptim0[ny]=m_ptSource[noff+1];
				ptim0[nh2-ny]=m_ptSource[noff+1];
			}
			tran.FFT(nh2,ptre0,ptim0,ptre1,ptim1);
			for(ny=0;ny<m_nHeight;ny++)
			{
				noff=ny*m_nWidth+nx;
				ptV[noff]=ptre1[ny]/nheight;
				ptV[noff+1]=ptim1[ny]/nheight;
			}
		}
		if(0!=nwidth%2)
		{
			for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=m_ptSource[ny*m_nWidth+nwidth];
			tran.EasyDCT(m_nHeight,ptIn,ptOut,false);
			for(ny=0;ny<m_nHeight;ny++) ptV[ny*m_nWidth+nwidth]=ptOut[ny];
		}
	}
	else for(nx=1;nx<=nwidth;nx++) for(ny=0;ny<m_nHeight;ny++) ptV[ny*m_nWidth+nx]=0;

	//Solve tridiagonal matrix equations
	for(ny=0;ny<m_nHeight;ny++)
	{	
		//Initialize tridiagonal matrix
		for(nx=2;nx<=nwidth;nx++) pta[nx]=1;
		for(nx=1;nx<=nwidth;nx++) ptb[nx]=-2-(float)pow(ny*tkappa,2.0);
		ptb[1]-=m_tBetaL/(m_tAlfaL*m_tdx-m_tBetaL);
		ptb[nwidth]+=m_tBetaH/(m_tAlfaH*m_tdx+m_tBetaH);
		for(nx=1;nx<nwidth;nx++) ptc[nx]=1;
		
		//Initialize right-hand side vector
		for(nx=1;nx<=nwidth;nx++) ptw[nx]=ptV[ny*m_nWidth+nx]*m_tdx*m_tdx;
		ptw[1]-=ptgamaL[ny]*m_tdx/(m_tAlfaL*m_tdx-m_tBetaL);
		ptw[nwidth]-=ptgamaH[ny]*m_tdx/(m_tAlfaH*m_tdx+m_tBetaH);
		
		//Invert tridiagonal matrix equation
		Tridiagonal(m_nWidth,pta,ptb,ptc,ptw,ptu);
		for(nx=1;nx<=nwidth;nx++) ptU[ny*m_nWidth+nx]=ptu[nx];
	}

	//Reconstruct solution via Fourier cosine transform
	for(nx=1;nx<nwidth;nx+=2)
	{
		ptre0[0]=ptU[nx];
		ptim0[0]=ptU[nx+1];
		noff=nheight*m_nWidth+nx;
		ptre0[nheight]=ptU[noff];
		ptim0[nheight]=ptU[noff+1];
		for(ny=1;ny<nheight;ny++)
		{
			noff=ny*m_nWidth+nx;
			ptre0[ny]=ptU[noff];
			ptre0[nh2-ny]=ptU[noff];
			ptim0[ny]=ptU[noff+1];
			ptim0[nh2-ny]=ptU[noff+1];
		}
		tran.FFT(nh2,ptre0,ptim0,ptre1,ptim1);
		for(ny=0;ny<m_nHeight;ny++)
		{
			noff=ny*m_nWidth+nx;
			m_ptSolver[noff]=ptre1[ny]/2;
			m_ptSolver[noff+1]=ptim1[ny]/2;
		}
	}
	if(0!=nwidth%2)
	{
		for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=ptU[ny*m_nWidth+nwidth];
		tran.EasyDCT(m_nHeight,ptIn,ptOut);
		for(ny=0;ny<m_nHeight;ny++) m_ptSolver[ny*m_nWidth+nwidth]=ptOut[ny];
	}

	//Calculate nx=0 and nx=m_nWidth-1 values
	for(ny=0;ny<m_nHeight;ny++)
	{
		m_ptSolver[ny*m_nWidth]=(m_ptBondL[ny]*m_tdx-m_tBetaL*m_ptSolver[ny*m_nWidth+1])/(m_tAlfaL*m_tdx-m_tBetaL);
		m_ptSolver[(ny+1)*m_nWidth-1]=(m_ptBondH[ny]*m_tdx+m_tBetaH*m_ptSolver[ny*m_nWidth+nwidth])/(m_tAlfaH*m_tdx+m_tBetaH);
	}

 	delete[] ptre0; delete[] ptim0;
	delete[] ptre1; delete[] ptim1;
	delete[] ptV; delete[] ptU;
	delete[] ptgamaL; delete[] ptgamaH;
	delete[] ptIn; delete[] ptOut;
	delete[] pta; delete[] ptb; delete[] ptc;
	delete[] ptw; delete[] ptu;

	return true;
}

CxPoisson2D::CxPoisson2D()
{
	m_nWidth=m_nHeight=0; m_tdx=m_tdy=0;
	m_txSuperior=m_tySuperior=1; m_txInferior=0;
	m_ptSource=m_ptSolver=NULL;
	m_tAlfaL=m_tBetaL=m_tAlfaH=m_tBetaH=0;
	m_ptBondL=m_ptBondH=NULL;
}

bool CxPoisson2D::Check1()
{
	if(0==m_nWidth || 0==m_nHeight) return false;
	if(0==m_tdx || 0==m_tdy) return false;
	if(NULL==m_ptSolver) return false;
	if((0!=m_tAlfaL || 0!=m_tBetaL) && NULL==m_ptBondL) return false;
	if((0!=m_tAlfaH || 0!=m_tBetaH) && NULL==m_ptBondH) return false;
	return true;
}

bool CxPoisson2D::Check2()
{
	if(fabs(m_tySuperior-m_tdy*(m_nHeight-1))>0.5*m_tdy) return false;
	if(fabs((m_txSuperior-m_txInferior)-m_tdx*(m_nWidth-1))>0.5*m_tdx) return false;
	return true;
}

//Function to invert tridiagonal matrix equation.Matrix is nlen by nlen.
//Left,centre and right diagonal elements is stored in pta,ptb,ptc.
//Right-hand side is stored in array ptw.Solution is written to array ptu.
//pta,ptb,ptc,ptw,ptu are of extent nlen+2 with redundant 0 and nlen+1 elements.
void CxPoisson2D::Tridiagonal(int ndim,double* pta,double* ptb,double* ptc,double* ptw,double*& ptu)
{
	int i,nlen=ndim-2;
	double *ptx=new double[nlen],*pty=new double[nlen];

	//Scan up diagonal from i=nlen to 1
	ptx[nlen-1]=-pta[nlen]/ptb[nlen];
	pty[nlen-1]=ptw[nlen]/ptb[nlen];
	for(i=nlen-2;i>0;i--)
	{
		ptx[i]=-pta[i+1]/(ptb[i+1]+ptc[i+1]*ptx[i+1]);
		pty[i]=(ptw[i+1]-ptc[i+1]*pty[i+1])/(ptb[i+1]+ptc[i+1]*ptx[i+1]);
	}
	ptx[0]=0; pty[0]=(ptw[1]-ptc[1]*pty[1])/(ptb[1]+ptc[1]*ptx[1]);

	//Scan down diagonal from i=1 to nlen
	ptu[1]=pty[0];
	for(i=1;i<nlen;i++) ptu[i+1]=ptx[i]*ptu[i]+pty[i];

	delete[] ptx; delete[] pty;
}

bool CxPoisson2D::Periodic()
{
	if(false==Check1()) return false;

	int nx,ny,nxid,nyid,noff,nsize=m_nWidth*m_nHeight;
	double txx,tyy,t2pi=8*atan(1);
	double *ptre=new double[nsize],*ptim=new double[nsize],*ptmp=new double[nsize];

	//Since the condition is periodic,the mean value of m_ptSource should be zero.
	fft2d.FFT2D(m_nWidth,m_nHeight,m_ptSource,NULL,ptre,ptim);
	for(ny=0;ny<m_nHeight;ny++)
	{
		nyid=ny>m_nHeight/2 ? ny-m_nHeight : ny;
		noff=ny*m_nWidth;
		for(nx=0;nx<m_nWidth;nx++)
		{
			if(0==nx && 0==ny) continue;
			nxid=nx>m_nWidth/2 ? nx-m_nWidth : nx;
			txx=2*(cos(t2pi*nxid/m_nWidth)-1)/pow(m_tdx,2.0);
			tyy=2*(cos(t2pi*nyid/m_nHeight)-1)/pow(m_tdy,2.0);
			ptre[noff+nx]/=txx+tyy; ptim[noff+nx]/=txx+tyy;
		}
	}
	fft2d.FFT2D(m_nWidth,m_nHeight,ptre,ptim,m_ptSolver,ptmp,false);
	delete[] ptre; delete[] ptim; delete[] ptmp;

	return true;
}

bool CxPoisson2D::EasyDirichlet()
{
	if(false==Check1() || false==Check2()) return false;

	int nx,ny,nwidth=m_nWidth-2,nheight=m_nHeight-1;
	double *ptV=new double[m_nWidth*m_nHeight],*ptU=new double[m_nWidth*m_nHeight];
	double *ptgamaL=new double[m_nHeight],*ptgamaH=new double[m_nHeight];
	double *ptIn=new double[m_nHeight],*ptOut=new double[m_nHeight];
	double *pta=new double[m_nWidth],*ptb=new double[m_nWidth],*ptc=new double[m_nWidth];
	double tpi=4*atan(1),*ptw=new double[m_nWidth],*ptu=new double[m_nWidth];
	double tkappa=tpi*m_tdx/m_tySuperior;

	//Inverse Fourier sine transform boundary conditions
	tran.EasyDST0(m_nHeight,m_ptBondL,ptgamaL,false);
	tran.EasyDST0(m_nHeight,m_ptBondH,ptgamaH,false);

	for(nx=1;nx<=nwidth;nx++)
	{
		if(NULL!=m_ptSource)//Inverse Fourier sine transform source term
		{
			for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=m_ptSource[ny*m_nWidth+nx];
			tran.EasyDST0(m_nHeight,ptIn,ptOut,false);
			for(ny=0;ny<m_nHeight;ny++) ptV[ny*m_nWidth+nx]=ptOut[ny];
		}
		else for(ny=0;ny<m_nHeight;ny++) ptV[ny*m_nWidth+nx]=0;
	}

	for(ny=1;ny<nheight;ny++)//Solve tridiagonal matrix equations
	{
		//Initialize tridiagonal matrix
		for(nx=2;nx<=nwidth;nx++) pta[nx]=1;
		for(nx=1;nx<=nwidth;nx++) ptb[nx]=-2-pow(ny*tkappa,2.0);
		ptb[1]-=m_tBetaL/(m_tAlfaL*m_tdx-m_tBetaL);
		ptb[nwidth]+=m_tBetaH/(m_tAlfaH*m_tdx+m_tBetaH);
		for(nx=1;nx<nwidth;nx++) ptc[nx]=1;

		//Initialize right-hand side vector
		for(nx=1;nx<=nwidth;nx++) ptw[nx]=ptV[ny*m_nWidth+nx]*m_tdx*m_tdx;
		ptw[1]-=ptgamaL[ny]*m_tdx/(m_tAlfaL*m_tdx-m_tBetaL);
		ptw[nwidth]-=ptgamaH[ny]*m_tdx/(m_tAlfaH*m_tdx+m_tBetaH);

		//Invert tridiagonal matrix equation
		Tridiagonal(m_nWidth,pta,ptb,ptc,ptw,ptu);
		for(nx=1;nx<=nwidth;nx++) ptU[ny*m_nWidth+nx]=ptu[nx];
	}

	for(nx=1;nx<=nwidth;nx++) ptU[nx]=ptU[nheight*m_nWidth+nx]=0;

	for(nx=1;nx<=nwidth;nx++)//Reconstruct solution via Fourier sine transform
	{
		for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=ptU[ny*m_nWidth+nx];
		tran.EasyDST0(m_nHeight,ptIn,ptOut);
		for(ny=0;ny<m_nHeight;ny++) m_ptSolver[ny*m_nWidth+nx]=ptOut[ny];
	}

	for(ny=0;ny<m_nHeight;ny++)//Calculate nx=0 and nx=m_nWidth-1 values
	{
		m_ptSolver[ny*m_nWidth]=(m_ptBondL[ny]*m_tdx-m_tBetaL*m_ptSolver[ny*m_nWidth+1])/(m_tAlfaL*m_tdx-m_tBetaL);
		m_ptSolver[(ny+1)*m_nWidth-1]=(m_ptBondH[ny]*m_tdx+m_tBetaH*m_ptSolver[ny*m_nWidth+nwidth])/(m_tAlfaH*m_tdx+m_tBetaH);
	}

	delete[] ptV; delete[] ptU;
	delete[] ptgamaL; delete[] ptgamaH;
	delete[] ptIn; delete[] ptOut;
	delete[] pta; delete[] ptb; delete[] ptc;
	delete[] ptw; delete[] ptu;

	return true;
}

bool CxPoisson2D::EasyNeumann()
{
	if(false==Check1() || false==Check2()) return false;

	int nx,ny,nwidth=m_nWidth-2,nheight=m_nHeight-1;
	double *ptV=new double[m_nWidth*m_nHeight],*ptU=new double[m_nWidth*m_nHeight];
	double *ptgamaL=new double[m_nHeight],*ptgamaH=new double[m_nHeight];
	double *ptIn=new double[m_nHeight],*ptOut=new double[m_nHeight];
	double *pta=new double[m_nWidth],*ptb=new double[m_nWidth],*ptc=new double[m_nWidth];
	double tpi=4*atan(1),*ptw=new double[m_nWidth],*ptu=new double[m_nWidth];
	double tkappa=tpi*m_tdx/m_tySuperior;

	//Inverse Fourier cosine transform boundary conditions
	tran.EasyDCT(m_nHeight,m_ptBondL,ptgamaL,false);
	tran.EasyDCT(m_nHeight,m_ptBondH,ptgamaH,false);

	for(nx=1;nx<=nwidth;nx++)
	{
		if(NULL!=m_ptSource)//Inverse Fourier cosine transform source term
		{
			for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=m_ptSource[ny*m_nWidth+nx];
			tran.EasyDCT(m_nHeight,ptIn,ptOut,false);
			for(ny=0;ny<m_nHeight;ny++) ptV[ny*m_nWidth+nx]=ptOut[ny];
		}
		else for(ny=0;ny<m_nHeight;ny++) ptV[ny*m_nWidth+nx]=0;
	}

	for(ny=0;ny<m_nHeight;ny++)//Solve tridiagonal matrix equations
	{	
		//Initialize tridiagonal matrix
		for(nx=2;nx<=nwidth;nx++) pta[nx]=1;
		for(nx=1;nx<=nwidth;nx++) ptb[nx]=-2-pow(ny*tkappa,2.0);
		ptb[1]-=m_tBetaL/(m_tAlfaL*m_tdx-m_tBetaL);
		ptb[nwidth]+=m_tBetaH/(m_tAlfaH*m_tdx+m_tBetaH);
		for(nx=1;nx<nwidth;nx++) ptc[nx]=1;
		
		//Initialize right-hand side vector
		for(nx=1;nx<=nwidth;nx++) ptw[nx]=ptV[ny*m_nWidth+nx]*m_tdx*m_tdx;
		ptw[1]-=ptgamaL[ny]*m_tdx/(m_tAlfaL*m_tdx-m_tBetaL);
		ptw[nwidth]-=ptgamaH[ny]*m_tdx/(m_tAlfaH*m_tdx+m_tBetaH);
		
		//Invert tridiagonal matrix equation
		Tridiagonal(m_nWidth,pta,ptb,ptc,ptw,ptu);
		for(nx=1;nx<=nwidth;nx++) ptU[ny*m_nWidth+nx]=ptu[nx];
	}

	for(nx=1;nx<=nwidth;nx++)//Reconstruct solution via Fourier cosine transform
	{
		for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=ptU[ny*m_nWidth+nx];
		tran.EasyDCT(m_nHeight,ptIn,ptOut);
		for(ny=0;ny<m_nHeight;ny++) m_ptSolver[ny*m_nWidth+nx]=ptOut[ny];
	}

	for(ny=0;ny<m_nHeight;ny++)//Calculate nx=0 and nx=m_nWidth-1 values
	{
		m_ptSolver[ny*m_nWidth]=(m_ptBondL[ny]*m_tdx-m_tBetaL*m_ptSolver[ny*m_nWidth+1])/(m_tAlfaL*m_tdx-m_tBetaL);
		m_ptSolver[(ny+1)*m_nWidth-1]=(m_ptBondH[ny]*m_tdx+m_tBetaH*m_ptSolver[ny*m_nWidth+nwidth])/(m_tAlfaH*m_tdx+m_tBetaH);
	}

	delete[] ptV; delete[] ptU;
	delete[] ptgamaL; delete[] ptgamaH;
	delete[] ptIn; delete[] ptOut;
	delete[] pta; delete[] ptb; delete[] ptc;
	delete[] ptw; delete[] ptu;

	return true;
}

bool CxPoisson2D::FastDirichlet()
{
	if(false==Check1() || false==Check2()) return false;

	int nx,ny,nwidth=m_nWidth-2,nheight=m_nHeight-1;
	double *ptV=new double[m_nWidth*m_nHeight],*ptU=new double[m_nWidth*m_nHeight];
	double *ptgamaL=new double[m_nHeight],*ptgamaH=new double[m_nHeight];
	double *ptIn=new double[m_nHeight],*ptOut=new double[m_nHeight];
	double *pta=new double[m_nWidth],*ptb=new double[m_nWidth],*ptc=new double[m_nWidth];
	double tpi=4*atan(1),*ptw=new double[m_nWidth],*ptu=new double[m_nWidth];
	double tkappa=tpi*m_tdx/m_tySuperior;

	int noff,nh2=2*nheight;
	double *ptre0=new double[nh2],*ptim0=new double[nh2];
	double *ptre1=new double[nh2],*ptim1=new double[nh2];

	//Inverse Fourier sine transform boundary conditions
	ptre0[0]=ptre0[nheight]=0;
	ptim0[0]=ptim0[nheight]=0;
	for(ny=1;ny<nheight;ny++)
	{
		ptre0[ny]=-m_ptBondL[ny];
		ptre0[nh2-ny]=m_ptBondL[ny];
		ptim0[ny]=-m_ptBondH[ny];
		ptim0[nh2-ny]=m_ptBondH[ny];
	}
	tran.FFT(nh2,ptre0,ptim0,ptre1,ptim1);
	for(ny=0;ny<m_nHeight;ny++)
	{
		ptgamaL[ny]=ptim1[ny]/nheight;
		ptgamaH[ny]=-ptre1[ny]/nheight;
	}

	if(NULL!=m_ptSource)//Inverse Fourier sine transform source term
	{
		for(nx=1;nx<nwidth;nx+=2)
		{
			ptre0[0]=ptre0[nheight]=0;
			ptim0[0]=ptim0[nheight]=0;
			for(ny=1;ny<nheight;ny++)
			{
				noff=ny*m_nWidth+nx;
				ptre0[ny]=-m_ptSource[noff];
				ptre0[nh2-ny]=m_ptSource[noff];
				ptim0[ny]=-m_ptSource[noff+1];
				ptim0[nh2-ny]=m_ptSource[noff+1];
			}
			tran.FFT(nh2,ptre0,ptim0,ptre1,ptim1);
			for(ny=0;ny<m_nHeight;ny++)
			{
				noff=ny*m_nWidth+nx;
				ptV[noff]=ptim1[ny]/nheight;
				ptV[noff+1]=-ptre1[ny]/nheight;
			}
		}
		if(0!=nwidth%2)
		{
			for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=m_ptSource[ny*m_nWidth+nwidth];
			tran.EasyDST0(m_nHeight,ptIn,ptOut,false);
			for(ny=0;ny<m_nHeight;ny++) ptV[ny*m_nWidth+nwidth]=ptOut[ny];
		}
	}
	else for(nx=1;nx<=nwidth;nx++) for(ny=0;ny<m_nHeight;ny++) ptV[ny*m_nWidth+nx]=0;

	//Solve tridiagonal matrix equations
	for(ny=1;ny<nheight;ny++)
	{
		//Initialize tridiagonal matrix
		for(nx=2;nx<=nwidth;nx++) pta[nx]=1;
		for(nx=1;nx<=nwidth;nx++) ptb[nx]=-2-pow(ny*tkappa,2.0);
		ptb[1]-=m_tBetaL/(m_tAlfaL*m_tdx-m_tBetaL);
		ptb[nwidth]+=m_tBetaH/(m_tAlfaH*m_tdx+m_tBetaH);
		for(nx=1;nx<nwidth;nx++) ptc[nx]=1;

		//Initialize right-hand side vector
		for(nx=1;nx<=nwidth;nx++) ptw[nx]=ptV[ny*m_nWidth+nx]*m_tdx*m_tdx;
		ptw[1]-=ptgamaL[ny]*m_tdx/(m_tAlfaL*m_tdx-m_tBetaL);
		ptw[nwidth]-=ptgamaH[ny]*m_tdx/(m_tAlfaH*m_tdx+m_tBetaH);

		//Invert tridiagonal matrix equation
		Tridiagonal(m_nWidth,pta,ptb,ptc,ptw,ptu);
		for(nx=1;nx<=nwidth;nx++) ptU[ny*m_nWidth+nx]=ptu[nx];
	}
	for(nx=1;nx<=nwidth;nx++) ptU[nx]=ptU[nheight*m_nWidth+nx]=0;

	//Reconstruct solution via Fourier sine transform
	for(nx=1;nx<nwidth;nx+=2)
	{
		ptre0[0]=ptre0[nheight]=0;
		ptim0[0]=ptim0[nheight]=0;
		for(ny=1;ny<nheight;ny++)
		{
			noff=ny*m_nWidth+nx;
			ptre0[ny]=-ptU[noff];
			ptre0[nh2-ny]=ptU[noff];
			ptim0[ny]=-ptU[noff+1];
			ptim0[nh2-ny]=ptU[noff+1];
		}
		tran.FFT(nh2,ptre0,ptim0,ptre1,ptim1);
		for(ny=0;ny<m_nHeight;ny++)
		{
			noff=ny*m_nWidth+nx;
			m_ptSolver[noff]=ptim1[ny]/2;
			m_ptSolver[noff+1]=-ptre1[ny]/2;
		}
	}
	if(0!=nwidth%2)
	{
		for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=ptU[ny*m_nWidth+nwidth];
		tran.EasyDST0(m_nHeight,ptIn,ptOut);
		for(ny=0;ny<m_nHeight;ny++) m_ptSolver[ny*m_nWidth+nwidth]=ptOut[ny];
	}

	//Calculate nx=0 and nx=m_nWidth-1 values
	for(ny=0;ny<m_nHeight;ny++)
	{
		m_ptSolver[ny*m_nWidth]=(m_ptBondL[ny]*m_tdx-m_tBetaL*m_ptSolver[ny*m_nWidth+1])/(m_tAlfaL*m_tdx-m_tBetaL);
		m_ptSolver[(ny+1)*m_nWidth-1]=(m_ptBondH[ny]*m_tdx+m_tBetaH*m_ptSolver[ny*m_nWidth+nwidth])/(m_tAlfaH*m_tdx+m_tBetaH);
	}

 	delete[] ptre0; delete[] ptim0;
	delete[] ptre1; delete[] ptim1;
	delete[] ptV; delete[] ptU;
	delete[] ptgamaL; delete[] ptgamaH;
	delete[] ptIn; delete[] ptOut;
	delete[] pta; delete[] ptb; delete[] ptc;
	delete[] ptw; delete[] ptu;

	return true;
}

bool CxPoisson2D::FastNeumann()
{
	if(false==Check1() || false==Check2()) return false;

	int nx,ny,nwidth=m_nWidth-2,nheight=m_nHeight-1;
	double *ptV=new double[m_nWidth*m_nHeight],*ptU=new double[m_nWidth*m_nHeight];
	double *ptgamaL=new double[m_nHeight],*ptgamaH=new double[m_nHeight];
	double *ptIn=new double[m_nHeight],*ptOut=new double[m_nHeight];
	double *pta=new double[m_nWidth],*ptb=new double[m_nWidth],*ptc=new double[m_nWidth];
	double tpi=4*atan(1),*ptw=new double[m_nWidth],*ptu=new double[m_nWidth];
	double tkappa=tpi*m_tdx/m_tySuperior;

	//Inverse Fourier cosine transform boundary conditions
	int noff,nh2=2*nheight;
	double *ptre0=new double[nh2],*ptim0=new double[nh2];
	double *ptre1=new double[nh2],*ptim1=new double[nh2];

	for(ny=0;ny<nheight;ny++)
	{
		ptre0[ny]=m_ptBondL[ny];
		ptre0[nh2-1-ny]=m_ptBondL[ny+1];
		ptim0[ny]=m_ptBondH[ny];
		ptim0[nh2-1-ny]=m_ptBondH[ny+1];
	}
	tran.FFT(nh2,ptre0,ptim0,ptre1,ptim1);
	for(ny=0;ny<m_nHeight;ny++)
	{
		ptgamaL[ny]=ptre1[ny]/nheight;
		ptgamaH[ny]=ptim1[ny]/nheight;
	}

	if(NULL!=m_ptSource)//Inverse Fourier cosine transform source term
	{
		for(nx=1;nx<nwidth;nx+=2)
		{
			ptre0[0]=m_ptSource[nx];
			ptim0[0]=m_ptSource[nx+1];
			noff=nheight*m_nWidth+nx;
			ptre0[nheight]=m_ptSource[noff];
			ptim0[nheight]=m_ptSource[noff+1];
			for(ny=1;ny<nheight;ny++)
			{
				noff=ny*m_nWidth+nx;
				ptre0[ny]=m_ptSource[noff];
				ptre0[nh2-ny]=m_ptSource[noff];
				ptim0[ny]=m_ptSource[noff+1];
				ptim0[nh2-ny]=m_ptSource[noff+1];
			}
			tran.FFT(nh2,ptre0,ptim0,ptre1,ptim1);
			for(ny=0;ny<m_nHeight;ny++)
			{
				noff=ny*m_nWidth+nx;
				ptV[noff]=ptre1[ny]/nheight;
				ptV[noff+1]=ptim1[ny]/nheight;
			}
		}
		if(0!=nwidth%2)
		{
			for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=m_ptSource[ny*m_nWidth+nwidth];
			tran.EasyDCT(m_nHeight,ptIn,ptOut,false);
			for(ny=0;ny<m_nHeight;ny++) ptV[ny*m_nWidth+nwidth]=ptOut[ny];
		}
	}
	else for(nx=1;nx<=nwidth;nx++) for(ny=0;ny<m_nHeight;ny++) ptV[ny*m_nWidth+nx]=0;

	//Solve tridiagonal matrix equations
	for(ny=0;ny<m_nHeight;ny++)
	{	
		//Initialize tridiagonal matrix
		for(nx=2;nx<=nwidth;nx++) pta[nx]=1;
		for(nx=1;nx<=nwidth;nx++) ptb[nx]=-2-pow(ny*tkappa,2.0);
		ptb[1]-=m_tBetaL/(m_tAlfaL*m_tdx-m_tBetaL);
		ptb[nwidth]+=m_tBetaH/(m_tAlfaH*m_tdx+m_tBetaH);
		for(nx=1;nx<nwidth;nx++) ptc[nx]=1;
		
		//Initialize right-hand side vector
		for(nx=1;nx<=nwidth;nx++) ptw[nx]=ptV[ny*m_nWidth+nx]*m_tdx*m_tdx;
		ptw[1]-=ptgamaL[ny]*m_tdx/(m_tAlfaL*m_tdx-m_tBetaL);
		ptw[nwidth]-=ptgamaH[ny]*m_tdx/(m_tAlfaH*m_tdx+m_tBetaH);
		
		//Invert tridiagonal matrix equation
		Tridiagonal(m_nWidth,pta,ptb,ptc,ptw,ptu);
		for(nx=1;nx<=nwidth;nx++) ptU[ny*m_nWidth+nx]=ptu[nx];
	}

	//Reconstruct solution via Fourier cosine transform
	for(nx=1;nx<nwidth;nx+=2)
	{
		ptre0[0]=ptU[nx];
		ptim0[0]=ptU[nx+1];
		noff=nheight*m_nWidth+nx;
		ptre0[nheight]=ptU[noff];
		ptim0[nheight]=ptU[noff+1];
		for(ny=1;ny<nheight;ny++)
		{
			noff=ny*m_nWidth+nx;
			ptre0[ny]=ptU[noff];
			ptre0[nh2-ny]=ptU[noff];
			ptim0[ny]=ptU[noff+1];
			ptim0[nh2-ny]=ptU[noff+1];
		}
		tran.FFT(nh2,ptre0,ptim0,ptre1,ptim1);
		for(ny=0;ny<m_nHeight;ny++)
		{
			noff=ny*m_nWidth+nx;
			m_ptSolver[noff]=ptre1[ny]/2;
			m_ptSolver[noff+1]=ptim1[ny]/2;
		}
	}
	if(0!=nwidth%2)
	{
		for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=ptU[ny*m_nWidth+nwidth];
		tran.EasyDCT(m_nHeight,ptIn,ptOut);
		for(ny=0;ny<m_nHeight;ny++) m_ptSolver[ny*m_nWidth+nwidth]=ptOut[ny];
	}

	//Calculate nx=0 and nx=m_nWidth-1 values
	for(ny=0;ny<m_nHeight;ny++)
	{
		m_ptSolver[ny*m_nWidth]=(m_ptBondL[ny]*m_tdx-m_tBetaL*m_ptSolver[ny*m_nWidth+1])/(m_tAlfaL*m_tdx-m_tBetaL);
		m_ptSolver[(ny+1)*m_nWidth-1]=(m_ptBondH[ny]*m_tdx+m_tBetaH*m_ptSolver[ny*m_nWidth+nwidth])/(m_tAlfaH*m_tdx+m_tBetaH);
	}

 	delete[] ptre0; delete[] ptim0;
	delete[] ptre1; delete[] ptim1;
	delete[] ptV; delete[] ptU;
	delete[] ptgamaL; delete[] ptgamaH;
	delete[] ptIn; delete[] ptOut;
	delete[] pta; delete[] ptb; delete[] ptc;
	delete[] ptw; delete[] ptu;

	return true;
}
