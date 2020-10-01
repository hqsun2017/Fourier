//Filter.h : main header file

#include <math.h>

#if !defined(_FILTER_H_)
#define _FILTER_H_

double Distance(const int nx,const int ny,const int nwid,const int nhit);

double Buttorworth1(double fd,const double fd0,const double fn);

double Buttorworth2(double fd,const double fd0,const double fn);

double Exponential1(double fd,const double fd0,const double fn);

double Exponential2(double fd,const double fd0,const double fn);

double Trapezoid(double fd,const double fd0,const double fd1);

#endif //!defined(_FILTER_H_)
