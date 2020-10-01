#include "Color.h"

CColor::CColor()
{
	m_bGamma=false; m_fGamma=0.6;
	m_Colors[0]=RGB(0,0,0);
	m_Colors[1]=RGB(0,0,255);
	m_Colors[2]=RGB(102,102,102);
	m_Colors[3]=RGB(0,255,0);
	m_Colors[4]=RGB(153,153,153);
	m_Colors[5]=RGB(255,0,0);
	m_Colors[6]=RGB(204,204,204);
	m_Colors[7]=RGB(0,255,255);
	m_Colors[8]=RGB(255,0,255);
	m_Colors[9]=RGB(255,255,0);
	m_Colors[10]=RGB(255,255,255);
}

COLORREF CColor::CalcColor(double ftmp,bool bcolor)
{
	//if(TOTALCOLORS<2) return;
	COLORREF color1=m_Colors[0],color2=m_Colors[TOTALCOLORS-1];
	if(m_fMin==m_fMax) return (bcolor ? color2 : RGB(255,255,255));
	if(ftmp>=m_fMax) return (bcolor ? color2 : RGB(255,255,255));
	if(ftmp<=m_fMin) return (bcolor ? color1 : RGB(0,0,0));

	if(false==bcolor)
	{
		int ngray=int(255.0*(ftmp-m_fMin)/(m_fMax-m_fMin));
		if(ngray>255) ngray=255;
		if(ngray<0) ngray=0;
		return RGB(ngray,ngray,ngray);
	}

	double ratio,offset=(m_fMax-m_fMin)/(TOTALCOLORS-1);
	int order=int((ftmp-m_fMin)/offset);
	if(order>=TOTALCOLORS-1) return color2;
	if(order<0) return color1;

	color1=m_Colors[order],color2=m_Colors[order+1];
	ratio=(ftmp-m_fMin)/offset-order;//Gradient color ratio

	//Gradient color
	unsigned char red,green,blue;
	red=(unsigned char)(GetRValue(color2)*ratio+GetRValue(color1)*(1-ratio));
	green=(unsigned char)(GetGValue(color2)*ratio+GetGValue(color1)*(1-ratio));
	blue=(unsigned char)(GetBValue(color2)*ratio+GetBValue(color1)*(1-ratio));
	if(true==m_bGamma)
	{
		red=(unsigned char)(pow((double)red/255.0,m_fGamma)*255);
		green=(unsigned char)(pow((double)green/255.0,m_fGamma)*255);
		blue=(unsigned char)(pow((double)blue/255.0,m_fGamma)*255);
	}

	return RGB(red,green,blue);
}
