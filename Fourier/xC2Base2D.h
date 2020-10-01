//xC2Base2D.h : main header file

#if !defined(_XC2BASE2D_H_)
#define _XC2BASE2D_H_

class CxC2Base2D_s
{
public:
	CxC2Base2D_s() {;}
	~CxC2Base2D_s() {;}

public:
	//循环运算
	virtual void CircleProcess(int nw,int nh,float* ptx,float* pty,float* ptc)=0;
	//线性运算
	void LinearProcess(int nw1,int nh1,float* ptx,int nw2,int nh2,float* pty,float* ptc);
};

class CxC2Base2D
{
public:
	CxC2Base2D() {;}
	~CxC2Base2D() {;}

public:
	//循环运算
	virtual void CircleProcess(int nw,int nh,double* ptx,double* pty,double* ptc)=0;
	//线性运算
	void LinearProcess(int nw1,int nh1,double* ptx,int nw2,int nh2,double* pty,double* ptc);
};

#endif //!defined(_XC2BASE2D_H_)
