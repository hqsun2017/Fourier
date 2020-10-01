//xC2Base1D.h : main header file

#if !defined(_XC2BASE1D_H_)
#define _XC2BASE1D_H_

class CxC2Base1D_s
{
public:
	CxC2Base1D_s() {;}
	~CxC2Base1D_s() {;}

public:
	//循环运算
	virtual void CircleProcess(int nSize,float* ptx,float* pty,float* ptc)=0;
	//线性运算
	void LinearProcess(int nSizx,float* ptx,int nSizy,float* pty,float* ptc);
};

class CxC2Base1D
{
public:
	CxC2Base1D() {;}
	~CxC2Base1D() {;}

public:
	//循环运算
	virtual void CircleProcess(int nSize,double* ptx,double* pty,double* ptc)=0;
	//线性运算
	void LinearProcess(int nSizx,double* ptx,int nSizy,double* pty,double* ptc);
};

#endif //!defined(_XC2BASE1D_H_)
