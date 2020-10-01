//xC2Base1D.h : main header file

#if !defined(_XC2BASE1D_H_)
#define _XC2BASE1D_H_

class CxC2Base1D_s
{
public:
	CxC2Base1D_s() {;}
	~CxC2Base1D_s() {;}

public:
	//ѭ������
	virtual void CircleProcess(int nSize,float* ptx,float* pty,float* ptc)=0;
	//��������
	void LinearProcess(int nSizx,float* ptx,int nSizy,float* pty,float* ptc);
};

class CxC2Base1D
{
public:
	CxC2Base1D() {;}
	~CxC2Base1D() {;}

public:
	//ѭ������
	virtual void CircleProcess(int nSize,double* ptx,double* pty,double* ptc)=0;
	//��������
	void LinearProcess(int nSizx,double* ptx,int nSizy,double* pty,double* ptc);
};

#endif //!defined(_XC2BASE1D_H_)
