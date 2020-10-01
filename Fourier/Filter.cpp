#include "Filter.h"

//To calculate the distance from (nx,ny) to (0,0) in the image of nwid*nhit
double Distance(const int nx,const int ny,const int nwid,const int nhit)
{
	double fx=nwid,fy=nhit;
	fx=fx/2-nx; fy=fy/2-ny;
	return sqrt(fx*fx+fy*fy);
}

double Buttorworth1(double fd,const double fd0,const double fn)
{
	if(0==fd0) return 0;
	return 1.0/(1.0+pow(fd/fd0,2.0*fn));
}

double Buttorworth2(double fd,const double fd0,const double fn)
{
	if(0==fd0) return 0;
	return 1.0/(1.0+(sqrt(2.0)-1)*pow(fd/fd0,2.0*fn));
}

double Exponential1(double fd,const double fd0,const double fn)
{
	if(0==fd0) return 0;
	return exp(-pow(fd/fd0,2.0*fn));
}

double Exponential2(double fd,const double fd0,const double fn)
{
	if(0==fd0) return 0;
	return exp(-log(2)/2.0*pow(fd/fd0,2.0*fn));
}

double Trapezoid(double fd,const double fd0,const double fd1)
{
	if(fd<fd0) return 1;
	else if(fd>fd1) return 0;
	else
	{
		if(fd1==fd0) return 0;
		else return (fd1-fd)/(fd1-fd0);
	}
}
