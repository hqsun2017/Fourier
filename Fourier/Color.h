//Color.h : main header file

#include <math.h>
#include <Windows.h>

#if !defined(_COLORS_H_)
#define _COLORS_H_

#define TOTALCOLORS 11

class CColor
{
public:
	CColor();
	~CColor() {;}

private:
	bool m_bGamma;
	double m_fGamma,m_fMax,m_fMin;
	COLORREF m_Colors[TOTALCOLORS];

public:
	bool& GammaState() {return m_bGamma;}
	double& GammaValue() {return m_fGamma;}
	double& MaxValue() {return m_fMax;}
	double& MinValue() {return m_fMin;}
	COLORREF CalcColor(double ftmp,bool bcolor=true);
};

#endif //!defined(_COLORS_H_)
