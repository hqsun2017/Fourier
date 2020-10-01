#include "xDiffusion2D.h"

CxDiffusion2D_s::CxDiffusion2D_s()
{
	m_nWidth=m_nHeight=0; m_tdt=m_tdx=m_tdy=0;
	m_txSuperior=m_tySuperior=1; m_txInferior=0;
	m_ptSource=m_ptSolver=NULL; m_ptBondL=m_ptBondH=NULL;
	m_tAlfaL=m_tBetaL=m_tAlfaH=m_tBetaH=0;
}

bool CxDiffusion2D_s::Check1()
{
	if(0==m_nWidth || 0==m_nHeight) return false;
	if(0==m_tD || 0==m_tdt || 0==m_tdx || 0==m_tdy) return false;
	if(NULL==m_ptSolver) return false;
	if((0!=m_tAlfaL || 0!=m_tBetaL) && NULL==m_ptBondL) return false;
	if((0!=m_tAlfaH || 0!=m_tBetaH) && NULL==m_ptBondH) return false;
	return true;
}

bool CxDiffusion2D_s::Check2()
{
	if(fabs(m_tySuperior-m_tdy*(m_nHeight-1))>0.5*m_tdy) return false;
	if(fabs((m_txSuperior-m_txInferior)-m_tdx*(m_nWidth-1))>0.5*m_tdx) return false;
	return true;
}

//Function to invert tridiagonal matrix equation.Matrix is nlen by nlen.
//Left,centre and right diagonal elements is stored in pta,ptb,ptc.
//Right-hand side is stored in array ptw.Solution is written to array ptu.
//pta,ptb,ptc,ptw,ptu are of extent nlen+2 with redundant 0 and nlen+1 elements.
void CxDiffusion2D_s::Tridiagonal(int ndim,float* pta,float* ptb,float* ptc,float* ptw,float*& ptu)
{
	int i,nlen=ndim-2;
	float *ptx=new float[nlen],*pty=new float[nlen];

	//Scan up diagonal from i=nlen to 1
	ptx[nlen-1]=-pta[nlen]/ptb[nlen]; pty[nlen-1]=ptw[nlen]/ptb[nlen];
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

bool CxDiffusion2D_s::Dirichlet()
{
	if(false==Check1() || false==Check2()) return false;

	int nx,ny,nwidth=m_nWidth-2,nheight=m_nHeight-1;
	float *ptF,*ptV=new float[m_nWidth*m_nHeight],*ptT=new float[m_nWidth*m_nHeight];
	float *ptgamaL=new float[m_nHeight],*ptgamaH=new float[m_nHeight];
	float tpi=4*(float)atan(1),*ptIn=new float[m_nHeight],*ptOut=new float[m_nHeight];
	float *pta=new float[m_nWidth],*ptb=new float[m_nWidth],*ptc=new float[m_nWidth];
	float *ptw=new float[m_nWidth],*ptu=new float[m_nWidth];
	float tC=0.5f*m_tD*m_tdt/(m_tdx*m_tdx),tkappa=tpi*m_tdx/m_tySuperior;

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

	//Inverse Fourier sine transform float	
	for(nx=0;nx<m_nWidth/2;nx++)
	{
		ptre0[0]=ptre0[nheight]=0;
		ptim0[0]=ptim0[nheight]=0;
		for(ny=1;ny<nheight;ny++)
		{
			noff=ny*m_nWidth+2*nx;
			ptre0[ny]=-m_ptSolver[noff];
			ptre0[nh2-ny]=m_ptSolver[noff];
			ptim0[ny]=-m_ptSolver[noff+1];
			ptim0[nh2-ny]=m_ptSolver[noff+1];
		}
		tran.FFT(nh2,ptre0,ptim0,ptre1,ptim1);
		for(ny=0;ny<m_nHeight;ny++)
		{
			noff=ny*m_nWidth+2*nx;
			ptT[noff]=ptim1[ny]/nheight;
			ptT[noff+1]=-ptre1[ny]/nheight;
		}
	}
	if(0!=m_nWidth%2)
	{
		for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=m_ptSolver[(ny+1)*m_nWidth-1];
		tran.EasyDST0(m_nHeight,ptIn,ptOut,false);
		for(ny=0;ny<m_nHeight;ny++) ptT[(ny+1)*m_nWidth-1]=ptOut[ny];
	}

	//Inverse Fourier sine transform source term
	if(NULL!=m_ptSource)
	{
		ptF=new float[m_nWidth*m_nHeight];
		for(nx=0;nx<m_nWidth/2;nx++)
		{
			ptre0[0]=ptre0[nheight]=0;
			ptim0[0]=ptim0[nheight]=0;
			for(ny=1;ny<nheight;ny++)
			{
				noff=ny*m_nWidth+2*nx;
				ptre0[ny]=-m_ptSource[noff];
				ptre0[nh2-ny]=m_ptSource[noff];
				ptim0[ny]=-m_ptSource[noff+1];
				ptim0[nh2-ny]=m_ptSource[noff+1];
			}
			tran.FFT(nh2,ptre0,ptim0,ptre1,ptim1);
			for(ny=0;ny<m_nHeight;ny++)
			{
				noff=ny*m_nWidth+2*nx;
				ptF[noff]=ptim1[ny]/nheight;
				ptF[noff+1]=-ptre1[ny]/nheight;
			}
		}
		if(0!=m_nWidth%2)
		{
			for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=m_ptSource[(ny+1)*m_nWidth-1];
			tran.EasyDST0(m_nHeight,ptIn,ptOut,false);
			for(ny=0;ny<m_nHeight;ny++) ptF[(ny+1)*m_nWidth-1]=ptOut[ny];
		}
	}

	//Construct source term
	for(nx=1;nx<=nwidth;nx++) for(ny=0;ny<m_nHeight;ny++)
	{
		noff=ny*m_nWidth+nx;
		ptV[noff]=tC*ptT[noff-1]+(1-tC*(2+(float)pow(ny*tkappa,2.0)))*ptT[noff]+tC*ptT[noff+1];
	}

	//Solve tridiagonal matrix equations
	for(ny=1;ny<nheight;ny++)
	{
		//Initialize tridiagonal matrix
		for(nx=2;nx<=nwidth;nx++) pta[nx]=-tC;
		for(nx=1;nx<=nwidth;nx++) ptb[nx]=1+tC*(2+(float)pow(ny*tkappa,2.0));
		ptb[1]+=tC*m_tBetaL/(m_tAlfaL*m_tdx-m_tBetaL);
		ptb[nwidth]-=tC*m_tBetaH/(m_tAlfaH*m_tdx+m_tBetaH);
		for(nx=1;nx<nwidth;nx++) ptc[nx]=-tC;

		//Initialize right-hand side vector
		for(nx=1;nx<=nwidth;nx++)
		{
			noff=ny*m_nWidth+nx;
			ptw[nx]=ptV[noff];
			if(NULL!=m_ptSource) ptw[nx]+=m_tdt*ptF[noff];
		}
		ptw[1]+=tC*ptgamaL[ny]*m_tdx/(m_tAlfaL*m_tdx-m_tBetaL);
		ptw[nwidth]+=tC*ptgamaH[ny]*m_tdx/(m_tAlfaH*m_tdx+m_tBetaH);

		//Invert tridiagonal matrix equation
		Tridiagonal(m_nWidth,pta,ptb,ptc,ptw,ptu);
		for(nx=1;nx<=nwidth;nx++) ptT[ny*m_nWidth+nx]=ptu[nx];
	}
	for(nx=1;nx<=nwidth;nx++) ptT[nx]=ptT[nheight*m_nWidth+nx]=0;

	//Reconstruct solution via Fourier sine transform
	for(nx=1;nx<nwidth;nx+=2)
	{
		ptre0[0]=ptre0[nheight]=0;
		ptim0[0]=ptim0[nheight]=0;
		for(ny=1;ny<nheight;ny++)
		{
			noff=ny*m_nWidth+nx;
			ptre0[ny]=-ptT[noff];
			ptre0[nh2-ny]=ptT[noff];
			ptim0[ny]=-ptT[noff+1];
			ptim0[nh2-ny]=ptT[noff+1];
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
		for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=ptT[ny*m_nWidth+nwidth];
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
	delete[] ptV; delete[] ptT;
	if(NULL!=m_ptSource) delete[] ptF;
	delete[] ptgamaL; delete[] ptgamaH;
	delete[] ptIn; delete[] ptOut;
	delete[] pta; delete[] ptb; delete[] ptc;
	delete[] ptw; delete[] ptu;

	return true;
}

bool CxDiffusion2D_s::Neumann()
{
	if(false==Check1() || false==Check2()) return false;

	int nx,ny,nwidth=m_nWidth-2,nheight=m_nHeight-1;
	float *ptF,*ptV=new float[m_nWidth*m_nHeight],*ptT=new float[m_nWidth*m_nHeight];
	float *ptgamaL=new float[m_nHeight],*ptgamaH=new float[m_nHeight];
	float tpi=4*(float)atan(1),*ptIn=new float[m_nHeight],*ptOut=new float[m_nHeight];
	float *pta=new float[m_nWidth],*ptb=new float[m_nWidth],*ptc=new float[m_nWidth];
	float *ptw=new float[m_nWidth],*ptu=new float[m_nWidth];
	float tC=0.5f*m_tD*m_tdt/(m_tdx*m_tdx),tkappa=tpi*m_tdx/m_tySuperior;

	int noff,nh2=2*nheight;
	float *ptre0=new float[nh2],*ptim0=new float[nh2];
	float *ptre1=new float[nh2],*ptim1=new float[nh2];

	//Inverse Fourier cosine transform boundary conditions
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

	//Inverse Fourier cosine transform float	
	for(nx=0;nx<m_nWidth/2;nx++)
	{
		ptre0[0]=m_ptSolver[2*nx];
		ptim0[0]=m_ptSolver[2*nx+1];
		noff=nheight*m_nWidth+2*nx;
		ptre0[nheight]=m_ptSolver[noff];
		ptim0[nheight]=m_ptSolver[noff+1];
		for(ny=1;ny<nheight;ny++)
		{
			noff=ny*m_nWidth+2*nx;
			ptre0[ny]=m_ptSolver[noff];
			ptre0[nh2-ny]=m_ptSolver[noff];
			ptim0[ny]=m_ptSolver[noff+1];
			ptim0[nh2-ny]=m_ptSolver[noff+1];
		}
		tran.FFT(nh2,ptre0,ptim0,ptre1,ptim1);
		for(ny=0;ny<m_nHeight;ny++)
		{
			noff=ny*m_nWidth+2*nx;
			ptT[noff]=ptre1[ny]/nheight;
			ptT[noff+1]=ptim1[ny]/nheight;
		}
	}
	if(0!=m_nWidth%2)
	{
		for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=m_ptSolver[(ny+1)*m_nWidth-1];
		tran.EasyDCT(m_nHeight,ptIn,ptOut,false);
		for(ny=0;ny<m_nHeight;ny++) ptT[(ny+1)*m_nWidth-1]=ptOut[ny];
	}

	//Inverse Fourier sine transform source term
	if(NULL!=m_ptSource)
	{
		ptF=new float[m_nWidth*m_nHeight];
		for(nx=0;nx<m_nWidth/2;nx++)
		{
			ptre0[0]=m_ptSource[2*nx];
			ptim0[0]=m_ptSource[2*nx+1];
			noff=nheight*m_nWidth+2*nx;
			ptre0[nheight]=m_ptSource[noff];
			ptim0[nheight]=m_ptSource[noff+1];
			for(ny=1;ny<nheight;ny++)
			{
				noff=ny*m_nWidth+2*nx;
				ptre0[ny]=m_ptSource[noff];
				ptre0[nh2-ny]=m_ptSource[noff];
				ptim0[ny]=m_ptSource[noff+1];
				ptim0[nh2-ny]=m_ptSource[noff+1];
			}
			tran.FFT(nh2,ptre0,ptim0,ptre1,ptim1);
			for(ny=0;ny<m_nHeight;ny++)
			{
				noff=ny*m_nWidth+2*nx;
				ptF[noff]=ptre1[ny]/nheight;
				ptF[noff+1]=ptim1[ny]/nheight;
			}
		}
		if(0!=m_nWidth%2)
		{
			for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=m_ptSource[(ny+1)*m_nWidth-1];
			tran.EasyDCT(m_nHeight,ptIn,ptOut,false);
			for(ny=0;ny<m_nHeight;ny++) ptF[(ny+1)*m_nWidth-1]=ptOut[ny];
		}
	}

	//Construct source term
	for(nx=1;nx<=nwidth;nx++) for(ny=0;ny<m_nHeight;ny++)
	{
		noff=ny*m_nWidth+nx;
		ptV[noff]=tC*ptT[noff-1]+(1-tC*(2+(float)pow(ny*tkappa,2)))*ptT[noff]+tC*ptT[noff+1];
	}

	//Solve tridiagonal matrix equations
	for(ny=0;ny<m_nHeight;ny++)
	{
		//Initialize tridiagonal matrix
		for(nx=2;nx<=nwidth;nx++) pta[nx]=-tC;
		for(nx=1;nx<=nwidth;nx++) ptb[nx]=1+tC*(2+(float)pow(ny*tkappa,2));
		ptb[1]+=tC*m_tBetaL/(m_tAlfaL*m_tdx-m_tBetaL);
		ptb[nwidth]-=tC*m_tBetaH/(m_tAlfaH*m_tdx+m_tBetaH);
		for(nx=1;nx<nwidth;nx++) ptc[nx]=-tC;

		//Initialize right-hand side vector
		for(nx=1;nx<=nwidth;nx++)
		{
			noff=ny*m_nWidth+nx;
			ptw[nx]=ptV[noff];
			if(NULL!=m_ptSource) ptw[nx]+=m_tdt*ptF[noff];
		}
		ptw[1]+=tC*ptgamaL[ny]*m_tdx/(m_tAlfaL*m_tdx-m_tBetaL);
		ptw[nwidth]+=tC*ptgamaH[ny]*m_tdx/(m_tAlfaH*m_tdx+m_tBetaH);

		//Invert tridiagonal matrix equation
		Tridiagonal(m_nWidth,pta,ptb,ptc,ptw,ptu);
		for(nx=1;nx<=nwidth;nx++) ptT[ny*m_nWidth+nx]=ptu[nx];
	}

	//Reconstruct solution via Fourier cosine transform
	for(nx=1;nx<nwidth;nx+=2)
	{
		ptre0[0]=ptT[nx];
		ptim0[0]=ptT[nx+1];
		noff=nheight*m_nWidth+nx;
		ptre0[nheight]=ptT[noff];
		ptim0[nheight]=ptT[noff+1];
		for(ny=1;ny<nheight;ny++)
		{
			noff=ny*m_nWidth+nx;
			ptre0[ny]=ptT[noff];
			ptre0[nh2-ny]=ptT[noff];
			ptim0[ny]=ptT[noff+1];
			ptim0[nh2-ny]=ptT[noff+1];
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
		for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=ptT[ny*m_nWidth+nwidth];
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
	delete[] ptV; delete[] ptT;
	if(NULL!=m_ptSource) delete[] ptF;
	delete[] ptgamaL; delete[] ptgamaH;
	delete[] ptIn; delete[] ptOut;
	delete[] pta; delete[] ptb; delete[] ptc;
	delete[] ptw; delete[] ptu;

	return true;
}

CxDiffusion2D::CxDiffusion2D()
{
	m_nWidth=m_nHeight=0; m_tdt=m_tdx=m_tdy=0;
	m_txSuperior=m_tySuperior=1; m_txInferior=0;
	m_ptSource=m_ptSolver=NULL; m_ptBondL=m_ptBondH=NULL;
	m_tAlfaL=m_tBetaL=m_tAlfaH=m_tBetaH=0;
}

bool CxDiffusion2D::Check1()
{
	if(0==m_nWidth || 0==m_nHeight) return false;
	if(0==m_tD || 0==m_tdt || 0==m_tdx || 0==m_tdy) return false;
	if(NULL==m_ptSolver) return false;
	if((0!=m_tAlfaL || 0!=m_tBetaL) && NULL==m_ptBondL) return false;
	if((0!=m_tAlfaH || 0!=m_tBetaH) && NULL==m_ptBondH) return false;
	return true;
}

bool CxDiffusion2D::Check2()
{
	if(fabs(m_tySuperior-m_tdy*(m_nHeight-1))>0.5*m_tdy) return false;
	if(fabs((m_txSuperior-m_txInferior)-m_tdx*(m_nWidth-1))>0.5*m_tdx) return false;
	return true;
}

//Function to invert tridiagonal matrix equation.Matrix is nlen by nlen.
//Left,centre and right diagonal elements is stored in pta,ptb,ptc.
//Right-hand side is stored in array ptw.Solution is written to array ptu.
//pta,ptb,ptc,ptw,ptu are of extent nlen+2 with redundant 0 and nlen+1 elements.
void CxDiffusion2D::Tridiagonal(int ndim,double* pta,double* ptb,double* ptc,double* ptw,double*& ptu)
{
	int i,nlen=ndim-2;
	double *ptx=new double[nlen],*pty=new double[nlen];

	//Scan up diagonal from i=nlen to 1
	ptx[nlen-1]=-pta[nlen]/ptb[nlen]; pty[nlen-1]=ptw[nlen]/ptb[nlen];
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

bool CxDiffusion2D::Dirichlet()
{
	if(false==Check1() || false==Check2()) return false;

	int nx,ny,nwidth=m_nWidth-2,nheight=m_nHeight-1;
	double *ptF,*ptV=new double[m_nWidth*m_nHeight],*ptT=new double[m_nWidth*m_nHeight];
	double *ptgamaL=new double[m_nHeight],*ptgamaH=new double[m_nHeight];
	double tpi=4*atan(1),*ptIn=new double[m_nHeight],*ptOut=new double[m_nHeight];
	double *pta=new double[m_nWidth],*ptb=new double[m_nWidth],*ptc=new double[m_nWidth];
	double *ptw=new double[m_nWidth],*ptu=new double[m_nWidth];
	double tC=0.5*m_tD*m_tdt/(m_tdx*m_tdx),tkappa=tpi*m_tdx/m_tySuperior;

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

	//Inverse Fourier sine transform double	
	for(nx=0;nx<m_nWidth/2;nx++)
	{
		ptre0[0]=ptre0[nheight]=0;
		ptim0[0]=ptim0[nheight]=0;
		for(ny=1;ny<nheight;ny++)
		{
			noff=ny*m_nWidth+2*nx;
			ptre0[ny]=-m_ptSolver[noff];
			ptre0[nh2-ny]=m_ptSolver[noff];
			ptim0[ny]=-m_ptSolver[noff+1];
			ptim0[nh2-ny]=m_ptSolver[noff+1];
		}
		tran.FFT(nh2,ptre0,ptim0,ptre1,ptim1);
		for(ny=0;ny<m_nHeight;ny++)
		{
			noff=ny*m_nWidth+2*nx;
			ptT[noff]=ptim1[ny]/nheight;
			ptT[noff+1]=-ptre1[ny]/nheight;
		}
	}
	if(0!=m_nWidth%2)
	{
		for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=m_ptSolver[(ny+1)*m_nWidth-1];
		tran.EasyDST0(m_nHeight,ptIn,ptOut,false);
		for(ny=0;ny<m_nHeight;ny++) ptT[(ny+1)*m_nWidth-1]=ptOut[ny];
	}

	//Inverse Fourier sine transform source term
	if(NULL!=m_ptSource)
	{
		ptF=new double[m_nWidth*m_nHeight];
		for(nx=0;nx<m_nWidth/2;nx++)
		{
			ptre0[0]=ptre0[nheight]=0;
			ptim0[0]=ptim0[nheight]=0;
			for(ny=1;ny<nheight;ny++)
			{
				noff=ny*m_nWidth+2*nx;
				ptre0[ny]=-m_ptSource[noff];
				ptre0[nh2-ny]=m_ptSource[noff];
				ptim0[ny]=-m_ptSource[noff+1];
				ptim0[nh2-ny]=m_ptSource[noff+1];
			}
			tran.FFT(nh2,ptre0,ptim0,ptre1,ptim1);
			for(ny=0;ny<m_nHeight;ny++)
			{
				noff=ny*m_nWidth+2*nx;
				ptF[noff]=ptim1[ny]/nheight;
				ptF[noff+1]=-ptre1[ny]/nheight;
			}
		}
		if(0!=m_nWidth%2)
		{
			for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=m_ptSource[(ny+1)*m_nWidth-1];
			tran.EasyDST0(m_nHeight,ptIn,ptOut,false);
			for(ny=0;ny<m_nHeight;ny++) ptF[(ny+1)*m_nWidth-1]=ptOut[ny];
		}
	}

	//Construct source term
	for(nx=1;nx<=nwidth;nx++) for(ny=0;ny<m_nHeight;ny++)
	{
		noff=ny*m_nWidth+nx;
		ptV[noff]=tC*ptT[noff-1]+(1-tC*(2+pow(ny*tkappa,2.0)))*ptT[noff]+tC*ptT[noff+1];
	}

	//Solve tridiagonal matrix equations
	for(ny=1;ny<nheight;ny++)
	{
		//Initialize tridiagonal matrix
		for(nx=2;nx<=nwidth;nx++) pta[nx]=-tC;
		for(nx=1;nx<=nwidth;nx++) ptb[nx]=1+tC*(2+pow(ny*tkappa,2.0));
		ptb[1]+=tC*m_tBetaL/(m_tAlfaL*m_tdx-m_tBetaL);
		ptb[nwidth]-=tC*m_tBetaH/(m_tAlfaH*m_tdx+m_tBetaH);
		for(nx=1;nx<nwidth;nx++) ptc[nx]=-tC;

		//Initialize right-hand side vector
		for(nx=1;nx<=nwidth;nx++)
		{
			noff=ny*m_nWidth+nx;
			ptw[nx]=ptV[noff];
			if(NULL!=m_ptSource) ptw[nx]+=m_tdt*ptF[noff];
		}
		ptw[1]+=tC*ptgamaL[ny]*m_tdx/(m_tAlfaL*m_tdx-m_tBetaL);
		ptw[nwidth]+=tC*ptgamaH[ny]*m_tdx/(m_tAlfaH*m_tdx+m_tBetaH);

		//Invert tridiagonal matrix equation
		Tridiagonal(m_nWidth,pta,ptb,ptc,ptw,ptu);
		for(nx=1;nx<=nwidth;nx++) ptT[ny*m_nWidth+nx]=ptu[nx];
	}
	for(nx=1;nx<=nwidth;nx++) ptT[nx]=ptT[nheight*m_nWidth+nx]=0;

	//Reconstruct solution via Fourier sine transform
	for(nx=1;nx<nwidth;nx+=2)
	{
		ptre0[0]=ptre0[nheight]=0;
		ptim0[0]=ptim0[nheight]=0;
		for(ny=1;ny<nheight;ny++)
		{
			noff=ny*m_nWidth+nx;
			ptre0[ny]=-ptT[noff];
			ptre0[nh2-ny]=ptT[noff];
			ptim0[ny]=-ptT[noff+1];
			ptim0[nh2-ny]=ptT[noff+1];
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
		for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=ptT[ny*m_nWidth+nwidth];
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
	delete[] ptV; delete[] ptT;
	if(NULL!=m_ptSource) delete[] ptF;
	delete[] ptgamaL; delete[] ptgamaH;
	delete[] ptIn; delete[] ptOut;
	delete[] pta; delete[] ptb; delete[] ptc;
	delete[] ptw; delete[] ptu;

	return true;
}

bool CxDiffusion2D::Neumann()
{
	if(false==Check1() || false==Check2()) return false;

	int nx,ny,nwidth=m_nWidth-2,nheight=m_nHeight-1;
	double *ptF,*ptV=new double[m_nWidth*m_nHeight],*ptT=new double[m_nWidth*m_nHeight];
	double *ptgamaL=new double[m_nHeight],*ptgamaH=new double[m_nHeight];
	double tpi=4*atan(1),*ptIn=new double[m_nHeight],*ptOut=new double[m_nHeight];
	double *pta=new double[m_nWidth],*ptb=new double[m_nWidth],*ptc=new double[m_nWidth];
	double *ptw=new double[m_nWidth],*ptu=new double[m_nWidth];
	double tC=0.5*m_tD*m_tdt/(m_tdx*m_tdx),tkappa=tpi*m_tdx/m_tySuperior;

	int noff,nh2=2*nheight;
	double *ptre0=new double[nh2],*ptim0=new double[nh2];
	double *ptre1=new double[nh2],*ptim1=new double[nh2];

	//Inverse Fourier cosine transform boundary conditions
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

	//Inverse Fourier cosine transform double	
	for(nx=0;nx<m_nWidth/2;nx++)
	{
		ptre0[0]=m_ptSolver[2*nx];
		ptim0[0]=m_ptSolver[2*nx+1];
		noff=nheight*m_nWidth+2*nx;
		ptre0[nheight]=m_ptSolver[noff];
		ptim0[nheight]=m_ptSolver[noff+1];
		for(ny=1;ny<nheight;ny++)
		{
			noff=ny*m_nWidth+2*nx;
			ptre0[ny]=m_ptSolver[noff];
			ptre0[nh2-ny]=m_ptSolver[noff];
			ptim0[ny]=m_ptSolver[noff+1];
			ptim0[nh2-ny]=m_ptSolver[noff+1];
		}
		tran.FFT(nh2,ptre0,ptim0,ptre1,ptim1);
		for(ny=0;ny<m_nHeight;ny++)
		{
			noff=ny*m_nWidth+2*nx;
			ptT[noff]=ptre1[ny]/nheight;
			ptT[noff+1]=ptim1[ny]/nheight;
		}
	}
	if(0!=m_nWidth%2)
	{
		for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=m_ptSolver[(ny+1)*m_nWidth-1];
		tran.EasyDCT(m_nHeight,ptIn,ptOut,false);
		for(ny=0;ny<m_nHeight;ny++) ptT[(ny+1)*m_nWidth-1]=ptOut[ny];
	}

	//Inverse Fourier sine transform source term
	if(NULL!=m_ptSource)
	{
		ptF=new double[m_nWidth*m_nHeight];
		for(nx=0;nx<m_nWidth/2;nx++)
		{
			ptre0[0]=m_ptSource[2*nx];
			ptim0[0]=m_ptSource[2*nx+1];
			noff=nheight*m_nWidth+2*nx;
			ptre0[nheight]=m_ptSource[noff];
			ptim0[nheight]=m_ptSource[noff+1];
			for(ny=1;ny<nheight;ny++)
			{
				noff=ny*m_nWidth+2*nx;
				ptre0[ny]=m_ptSource[noff];
				ptre0[nh2-ny]=m_ptSource[noff];
				ptim0[ny]=m_ptSource[noff+1];
				ptim0[nh2-ny]=m_ptSource[noff+1];
			}
			tran.FFT(nh2,ptre0,ptim0,ptre1,ptim1);
			for(ny=0;ny<m_nHeight;ny++)
			{
				noff=ny*m_nWidth+2*nx;
				ptF[noff]=ptre1[ny]/nheight;
				ptF[noff+1]=ptim1[ny]/nheight;
			}
		}
		if(0!=m_nWidth%2)
		{
			for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=m_ptSource[(ny+1)*m_nWidth-1];
			tran.EasyDCT(m_nHeight,ptIn,ptOut,false);
			for(ny=0;ny<m_nHeight;ny++) ptF[(ny+1)*m_nWidth-1]=ptOut[ny];
		}
	}

	//Construct source term
	for(nx=1;nx<=nwidth;nx++) for(ny=0;ny<m_nHeight;ny++)
	{
		noff=ny*m_nWidth+nx;
		ptV[noff]=tC*ptT[noff-1]+(1-tC*(2+pow(ny*tkappa,2)))*ptT[noff]+tC*ptT[noff+1];
	}

	//Solve tridiagonal matrix equations
	for(ny=0;ny<m_nHeight;ny++)
	{
		//Initialize tridiagonal matrix
		for(nx=2;nx<=nwidth;nx++) pta[nx]=-tC;
		for(nx=1;nx<=nwidth;nx++) ptb[nx]=1+tC*(2+pow(ny*tkappa,2));
		ptb[1]+=tC*m_tBetaL/(m_tAlfaL*m_tdx-m_tBetaL);
		ptb[nwidth]-=tC*m_tBetaH/(m_tAlfaH*m_tdx+m_tBetaH);
		for(nx=1;nx<nwidth;nx++) ptc[nx]=-tC;

		//Initialize right-hand side vector
		for(nx=1;nx<=nwidth;nx++)
		{
			noff=ny*m_nWidth+nx;
			ptw[nx]=ptV[noff];
			if(NULL!=m_ptSource) ptw[nx]+=m_tdt*ptF[noff];
		}
		ptw[1]+=tC*ptgamaL[ny]*m_tdx/(m_tAlfaL*m_tdx-m_tBetaL);
		ptw[nwidth]+=tC*ptgamaH[ny]*m_tdx/(m_tAlfaH*m_tdx+m_tBetaH);

		//Invert tridiagonal matrix equation
		Tridiagonal(m_nWidth,pta,ptb,ptc,ptw,ptu);
		for(nx=1;nx<=nwidth;nx++) ptT[ny*m_nWidth+nx]=ptu[nx];
	}

	//Reconstruct solution via Fourier cosine transform
	for(nx=1;nx<nwidth;nx+=2)
	{
		ptre0[0]=ptT[nx];
		ptim0[0]=ptT[nx+1];
		noff=nheight*m_nWidth+nx;
		ptre0[nheight]=ptT[noff];
		ptim0[nheight]=ptT[noff+1];
		for(ny=1;ny<nheight;ny++)
		{
			noff=ny*m_nWidth+nx;
			ptre0[ny]=ptT[noff];
			ptre0[nh2-ny]=ptT[noff];
			ptim0[ny]=ptT[noff+1];
			ptim0[nh2-ny]=ptT[noff+1];
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
		for(ny=0;ny<m_nHeight;ny++) ptIn[ny]=ptT[ny*m_nWidth+nwidth];
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
	delete[] ptV; delete[] ptT;
	if(NULL!=m_ptSource) delete[] ptF;
	delete[] ptgamaL; delete[] ptgamaH;
	delete[] ptIn; delete[] ptOut;
	delete[] pta; delete[] ptb; delete[] ptc;
	delete[] ptw; delete[] ptu;

	return true;
}
