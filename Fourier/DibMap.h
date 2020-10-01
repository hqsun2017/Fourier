#ifndef _DIBITMAP_H_
#define _DIBITMAP_H_

#include <afxwin.h>

#define WIDTH(x) (((x)+3)/4*4)

class CDIBitmap
{
// Construction
public:
	CDIBitmap();
	CDIBitmap(int width,int height,int bits=24);
	CDIBitmap(const CDIBitmap* ptrbmp);
	virtual ~CDIBitmap();

// Attributes
private:
	BITMAPINFO* m_pInfo;
	BYTE* m_pPixels;
	CPalette* m_pPal;
	bool m_bIsColTrue;
	int m_nPalEntries,m_nPixels;

// Operations
protected:
	void DestroyBitmap();
	bool CreateNewPalette();//auto. made by Load() and CreateFromBitmap()
	void ClearPalette();//destroy the palette associated with this image

	BITMAPINFO* GetHeaderPtr() const;
	RGBQUAD* GetColorTablePtr() const;
	BYTE* GetPixelPtr() const;
	CPalette* GetPalette() const {return m_pPal;}

	int GetPalEntries() const;
	int GetPalEntries(BITMAPINFOHEADER& infoHeader) const;

public:
	int GetWidth() const;
	int GetHeight() const;
	void Resize(int nw,int nh);

	WORD GetColorCount() const;
	DWORD GetBitsPerPixel() const;
	CDIBitmap& operator =(const CDIBitmap& myBmp);
	bool operator ==(const CDIBitmap& myBmp);
	bool operator !=(const CDIBitmap& myBmp);
	BYTE* GetPixelByte(int nx,int ny);
	BYTE GetPixelValue(int nx,int ny);
	void SetPixelByte(int nx,int ny,BYTE* pByte,int nlen=1);
	void ColorToGrey();
	void SortGreyLevel();

	//draw the bitmap at the specified location
	virtual void DrawDIB(CDC* pDC,int x=0,int y=0);
	//draw the bitmap and stretch/compress it to the desired size
	virtual void DrawDIB(CDC* pDC,int x,int y,int width,int height);
	//draw parts of the dib into a given area of the DC
	virtual int DrawDIB(CDC* pDC,CRect& rectDC,CRect& rectDIB);
	//load a bitmap from disk
	virtual bool Load(CFile* pFile);
	virtual bool Load(const CString&);
	//save the bitmap to disk
	virtual bool Save(CFile* pFile);
	virtual bool Save(const CString&);
};

#endif // _DIBITMAP_H_
