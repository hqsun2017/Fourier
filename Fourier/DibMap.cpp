#include "DibMap.h"

CDIBitmap::CDIBitmap()
{
	m_pInfo=NULL; m_pPixels=NULL;
	m_pPal=NULL; m_bIsColTrue=false;
}

CDIBitmap::CDIBitmap(int width,int height,int bits)
{
	m_pPal=NULL; m_bIsColTrue=false;
	WORD wColors;
	switch(bits)//Set the appropriate number of colors
	{
	case 1: wColors=2; break;
	case 4: wColors=16; break;
	case 8: wColors=256; break;
	default : wColors=0; m_bIsColTrue=true; break;
	}

	DWORD dwWidth;
	//dwWidth=((((width*bits)+31)&~31)>>3);
	if(bits>8) dwWidth=WIDTH(width*3);
	else dwWidth=WIDTH(width);

	m_nPixels=dwWidth*height;
	m_pPixels=new BYTE[m_nPixels];
	if(!m_pPixels) throw TEXT("不能分配数据存储单元\n");

	m_pInfo=(BITMAPINFO*)new BYTE[sizeof(BITMAPINFOHEADER)+wColors*sizeof(RGBQUAD)];
	if(!m_pInfo) throw TEXT("不能分配BITMAPINFO结构单元\n");

	//Populate BITMAPINFO header info
	m_pInfo->bmiHeader.biSize=sizeof(BITMAPINFOHEADER);
	m_pInfo->bmiHeader.biWidth=width;
	m_pInfo->bmiHeader.biHeight=height;
	m_pInfo->bmiHeader.biPlanes=1;//////////////////////////

	m_pInfo->bmiHeader.biBitCount=bits;
	m_pInfo->bmiHeader.biCompression=BI_RGB;
	m_pInfo->bmiHeader.biSizeImage=m_nPixels;
	m_pInfo->bmiHeader.biXPelsPerMeter=0;
	m_pInfo->bmiHeader.biYPelsPerMeter=0;
	m_pInfo->bmiHeader.biClrUsed=0;
	m_pInfo->bmiHeader.biClrImportant=0;
}

CDIBitmap::CDIBitmap(const CDIBitmap* ptrbmp)
{
	m_pInfo=NULL; m_pPixels=NULL; m_pPal=NULL;

	WORD wColors=ptrbmp->GetColorCount();
	int nInfo=sizeof(BITMAPINFOHEADER)+wColors*sizeof(RGBQUAD);
	m_pInfo=(BITMAPINFO*)new BYTE[nInfo];
	memcpy(m_pInfo,ptrbmp->GetHeaderPtr(),nInfo);
	m_nPalEntries=GetPalEntries();

	m_nPixels=(ptrbmp->GetHeight())*(ptrbmp->GetWidth())*(ptrbmp->GetBitsPerPixel()/8);
	m_pPixels=new BYTE[m_nPixels];
	memcpy(m_pPixels,ptrbmp->GetPixelPtr(),m_nPixels);

	if(GetBitsPerPixel()<=8) m_bIsColTrue=false;
	else m_bIsColTrue=true;
}

CDIBitmap::~CDIBitmap()
{
	DestroyBitmap();
}

void CDIBitmap::DestroyBitmap()
{
	if(m_pInfo!=NULL) {delete[] (BYTE*)m_pInfo; m_pInfo=NULL;}
	if(m_pPixels!=NULL) {delete[] m_pPixels; m_pPixels=NULL;}
	ClearPalette();
}

bool CDIBitmap::CreateNewPalette()
{
	ASSERT(m_pInfo);

	ClearPalette();
	//We only need a palette,if there are <=256 colors.otherwise we would bomb the memory.
	if(m_pInfo->bmiHeader.biBitCount<=8)
	{
		m_pPal=new CPalette;

		m_nPalEntries=GetPalEntries();
		int nPalette=sizeof(LOGPALETTE)+sizeof(PALETTEENTRY)*m_nPalEntries;
		//Since the LOGPALETTE structure is open-ended,you must
		//dynamically allocate it,rather than using one off the stack.
		LOGPALETTE* pPal=(LOGPALETTE*)new BYTE[nPalette];
		RGBQUAD* pColorTab=GetColorTablePtr();
		pPal->palVersion=0x300;
		pPal->palNumEntries=m_nPalEntries;
		//Roll through the color table,and add each color to the logical palette.
		for(int ndx=0;ndx<m_nPalEntries;ndx++)
		{
			pPal->palPalEntry[ndx].peRed=pColorTab[ndx].rgbRed;
			pPal->palPalEntry[ndx].peGreen=pColorTab[ndx].rgbGreen;
			pPal->palPalEntry[ndx].peBlue=pColorTab[ndx].rgbBlue;
			pPal->palPalEntry[ndx].peFlags=0;
		}
		m_pPal->CreatePalette(pPal);
		delete[] (BYTE*)pPal;
	}

	return m_pPal ? true : false;
}

void CDIBitmap::ClearPalette()
{
	if(NULL!=m_pPal) delete m_pPal; m_pPal=NULL;
}

BITMAPINFO* CDIBitmap::GetHeaderPtr() const
{
	ASSERT(m_pInfo);//ASSERT(m_pPixels);
	return m_pInfo;
}

RGBQUAD* CDIBitmap::GetColorTablePtr() const
{
	ASSERT(m_pInfo);//ASSERT(m_pPixels);
	RGBQUAD* pColorTable=NULL;
	if(m_pInfo!=NULL)
	{
		int cOffset=sizeof(BITMAPINFOHEADER);
		pColorTable=(RGBQUAD*)(((BYTE*)(m_pInfo))+cOffset);
	}
	return pColorTable;
}

BYTE* CDIBitmap::GetPixelPtr() const
{
	ASSERT(m_pInfo); ASSERT(m_pPixels);
	return m_pPixels;
}

int CDIBitmap::GetWidth() const
{
	ASSERT(m_pInfo);
	return m_pInfo->bmiHeader.biWidth;
}

int CDIBitmap::GetHeight() const
{
	ASSERT(m_pInfo);
	return m_pInfo->bmiHeader.biHeight;
}

WORD CDIBitmap::GetColorCount() const
{
	ASSERT(m_pInfo);

	switch(m_pInfo->bmiHeader.biBitCount)
	{
	case 1: return 2;
	case 4: return 16;
	case 8: return 256;
	default: return 0;
	}
}

DWORD CDIBitmap::GetBitsPerPixel() const
{
	ASSERT(m_pInfo);
	return m_pInfo->bmiHeader.biBitCount;
}

void CDIBitmap::Resize(int nw,int nh)
{
	delete[] m_pPixels; m_pPixels=NULL;
	delete m_pPal; m_pPal=NULL;

	m_pInfo->bmiHeader.biWidth=nw;
	m_pInfo->bmiHeader.biHeight=nh;
	m_nPixels=WIDTH(nw*(GetBitsPerPixel()/8))*nh;
	m_pInfo->bmiHeader.biSizeImage=m_nPixels;
	m_pPixels=new BYTE[m_nPixels];
}

int CDIBitmap::GetPalEntries() const
{
	ASSERT(m_pInfo);
	return GetPalEntries(*(BITMAPINFOHEADER*)m_pInfo);
}

int CDIBitmap::GetPalEntries(BITMAPINFOHEADER& infoHeader) const
{
	int nReturn;
	if(infoHeader.biClrUsed==0) nReturn=(1<<infoHeader.biBitCount);
	else nReturn=infoHeader.biClrUsed;

	return nReturn;
}

void CDIBitmap::DrawDIB(CDC* pDC,int x,int y)
{
	DrawDIB(pDC,x,y,GetWidth(),GetHeight());
}

//DrawDib uses StretchDIBits to display the bitmap.
void CDIBitmap::DrawDIB(CDC* pDC,int x,int y,int width,int height)
{
	ASSERT(pDC);
	HDC hdc=pDC->GetSafeHdc();

	CPalette* pOldPal=NULL;

	CreateNewPalette();//new code
	if(m_pPal)
	{
		pOldPal=pDC->SelectPalette(m_pPal,false);
		pDC->RealizePalette();
		//Make sure to use the stretching mode best for color pictures
		pDC->SetStretchBltMode(COLORONCOLOR);
	}

	if(m_pInfo) StretchDIBits(hdc,x,y,width,height,0,0,GetWidth(),
		GetHeight(),GetPixelPtr(),GetHeaderPtr(),DIB_RGB_COLORS,SRCCOPY);

	if(m_pPal) pDC->SelectPalette(pOldPal,false);
}

int CDIBitmap::DrawDIB(CDC* pDC,CRect& rectDC,CRect& rectDIB)
{
	ASSERT(pDC);
	HDC hdc=pDC->GetSafeHdc();

	CPalette* pOldPal=NULL;

	CreateNewPalette();//new code
	if(m_pPal)
	{
		pOldPal=pDC->SelectPalette(m_pPal,false);
		pDC->RealizePalette();
		//Make sure to use the stretching mode best for color pictures
		pDC->SetStretchBltMode(COLORONCOLOR);
	}

	int nRet=0;

	if(m_pInfo)	nRet=SetDIBitsToDevice(hdc,//device
		rectDC.left,rectDC.top,//DestX&DestY
		rectDC.Width(),rectDC.Height(),//DestWidth&DestHeight
		rectDIB.left,//SrcX
		GetHeight()-rectDIB.top-rectDIB.Height(),//SrcY
		0,GetHeight(),//StartScan&NumScans
		GetPixelPtr(),GetHeaderPtr(),//color data&header data
		DIB_RGB_COLORS);//color usage

	if(m_pPal) pDC->SelectPalette(pOldPal,false);

	return nRet;
}

bool CDIBitmap::Load(CFile* pFile)
{
	ASSERT(pFile);

	bool fReturn=true;
	try
	{
		if(m_pInfo!=NULL) delete[] (BYTE*)m_pInfo; m_pInfo=NULL;
		if(m_pPixels!=NULL) delete[] m_pPixels; m_pPixels=NULL;
		//DWORD dwStart=pFile->GetPosition();

		//Check to make sure we have a bitmap. The first is 'B' and 'M'.
		BITMAPFILEHEADER fileHeader;
		pFile->Read(&fileHeader,sizeof(fileHeader));
		if(fileHeader.bfType!=0x4D42) throw TEXT("错误:文件格式出错,不是一个位图文件\n");

		BITMAPINFOHEADER infoHeader;
		pFile->Read(&infoHeader,sizeof(infoHeader));
		if(infoHeader.biSize!=sizeof(infoHeader)) throw TEXT("错误:不支持OS2下的PM BMP格式\n");

		//Store the sizes of the DIB structures
		int cColorTable,cInfo;

		m_nPalEntries=GetPalEntries(infoHeader);
		if(m_nPalEntries<=256)
		{
			m_bIsColTrue=false;
			cColorTable=m_nPalEntries*sizeof(RGBQUAD);
		}
		else
		{
			m_bIsColTrue=true;
			cColorTable=0;
		}
		cInfo=sizeof(BITMAPINFOHEADER)+cColorTable;
		m_nPixels=fileHeader.bfSize-fileHeader.bfOffBits;

		//Allocate space for a new bitmap info header,and copy
		//the info header that was loaded from the file. Read the
		//the file and store the results in the color table.
		m_pInfo=(BITMAPINFO*)new BYTE[cInfo];
		memcpy(m_pInfo,&infoHeader,sizeof(BITMAPINFOHEADER));
		pFile->Read(((BYTE*)m_pInfo)+sizeof(BITMAPINFOHEADER),cColorTable);

		//Allocate space for the pixel area,and load the pixel info from the file.
		m_pPixels=new BYTE[m_nPixels];
		//pFile->Seek(dwStart+fileHeader.bfOffBits,CFile::begin);
		pFile->Read(m_pPixels,m_nPixels);
		//CreateNewPalette();

#ifdef _DEBUG
	}
	catch(TCHAR* psz)
	{
		TRACE(psz);
#else
	}
	catch(TCHAR*)
	{
#endif
		fReturn=false;
	}

	return fReturn;
}

bool CDIBitmap::Load(const CString& strFilename)
{
	CFile file;
	if(file.Open(strFilename,CFile::modeRead)) return Load(&file);
	return false;
}

//Does not open or close pFile. Assumes caller will do it.
bool CDIBitmap::Save(CFile* pFile)
{
	ASSERT_VALID(pFile); ASSERT(m_pInfo); ASSERT(m_pPixels);

	BITMAPFILEHEADER bmfHdr;

	DWORD dwPadWidth=WIDTH(GetWidth()*(GetBitsPerPixel()/8));

	bmfHdr.bfType=0x4D42;
	//initialize to BitmapInfo size
	DWORD dwImageSize=m_pInfo->bmiHeader.biSize;
	WORD wColors=GetColorCount();//Add in palette size
	WORD wPaletteSize=(WORD)(wColors*sizeof(RGBQUAD));
	dwImageSize+=wPaletteSize;

	//Add in size of actual bit array
	dwImageSize+=dwPadWidth*GetHeight();
	m_pInfo->bmiHeader.biSizeImage=0;
	bmfHdr.bfSize=dwImageSize+sizeof(BITMAPFILEHEADER);
	bmfHdr.bfReserved1=bmfHdr.bfReserved2=0;
	bmfHdr.bfOffBits=(DWORD)sizeof(BITMAPFILEHEADER)+m_pInfo->bmiHeader.biSize+wPaletteSize;
	pFile->Write(&bmfHdr,sizeof(BITMAPFILEHEADER));

	pFile->Write(m_pInfo,sizeof(BITMAPINFO)+(wColors-1)*sizeof(RGBQUAD));
	pFile->WriteHuge(m_pPixels,DWORD(dwPadWidth*GetHeight()));

	return true;
}

bool CDIBitmap::Save(const CString& strFileName)
{
	ASSERT(!strFileName.IsEmpty());

	CFile File;

	if(!File.Open(strFileName,CFile::modeCreate|CFile::modeWrite))
	{
		TRACE1("CDIBitmap::Save():打开 %s 写操作失败.\n",LPCSTR(strFileName));
		return false;
	}

	return Save(&File);
}

CDIBitmap& CDIBitmap::operator =(const CDIBitmap& myBmp)
{
	DestroyBitmap();

	WORD wColors=myBmp.GetColorCount();
	int nInfo=sizeof(BITMAPINFOHEADER)+wColors*sizeof(RGBQUAD);
	m_pInfo=(BITMAPINFO*)new BYTE[nInfo];
	memcpy(m_pInfo,myBmp.GetHeaderPtr(),nInfo);
	m_nPalEntries=GetPalEntries();

	m_nPixels=(myBmp.GetHeight())*WIDTH(myBmp.GetWidth()*(myBmp.GetBitsPerPixel()/8));
	m_pPixels=new BYTE[m_nPixels];
	memcpy(m_pPixels,myBmp.GetPixelPtr(),m_nPixels);

	if(GetBitsPerPixel()<=8) m_bIsColTrue=false;
	else m_bIsColTrue=true;

	return *this;
}

bool CDIBitmap::operator ==(const CDIBitmap& myBmp)
{
	if(m_pInfo==NULL) return false;
	if(m_pPixels==NULL) return false;
	//if(m_pPal==NULL) return false;

	if(myBmp.GetHeaderPtr()==NULL) return false;
	if(myBmp.GetPixelPtr()==NULL) return false;
	//if(myBmp.GetPalette()==NULL) return false;

	WORD wColors=myBmp.GetColorCount();
	if(wColors!=GetColorCount()) return false;
	int nInfo=sizeof(BITMAPINFOHEADER)+wColors*sizeof(RGBQUAD);
	if(memcmp((char*)m_pInfo,(char*)myBmp.GetHeaderPtr(),nInfo)!=0) return false;

	//if(memcmp(m_pPixels,myBmp.GetPixelPtr(),m_nPixels)!=0) return false;

	return true;
}

bool CDIBitmap::operator !=(const CDIBitmap& myBmp)
{
	if(*this==myBmp) return false;
	else return true;
}

BYTE* CDIBitmap::GetPixelByte(int nx,int ny)
{
	ASSERT(m_pPixels);

	BYTE* ptr;

	if((nx<0) || (nx>=GetWidth()) || (ny<0) || (ny>=GetHeight())) ptr=NULL;
	else ptr=m_pPixels+((GetHeight()-1-ny)*WIDTH(GetWidth()*(GetBitsPerPixel()/8))+nx*(GetBitsPerPixel()/8));

	return ptr;
}

BYTE CDIBitmap::GetPixelValue(int nx,int ny)
{
	ASSERT(m_pPixels);//ASSERT(m_pInfo);

	BYTE bt;

	if((nx<0) || (nx>=GetWidth()) || (ny<0) || (ny>=GetHeight())) bt=0;
	else bt=m_pPixels[((GetHeight()-1-ny)*WIDTH(GetWidth()*(GetBitsPerPixel()/8))+nx*(GetBitsPerPixel()/8))];

	return bt;
}

void CDIBitmap::SetPixelByte(int nx,int ny,BYTE* pByte,int nlen)
{
	ASSERT(m_pPixels);

	BYTE* ptr;

	if((nx<0) || (nx>=GetWidth()) || (ny<0) || (ny>=GetHeight())) return;
	else ptr=m_pPixels+((GetHeight()-1-ny)*WIDTH(GetWidth()*(GetBitsPerPixel()/8))+nx*(GetBitsPerPixel()/8));

	memcpy(ptr,pByte,nlen);
}

void CDIBitmap::ColorToGrey()
{
	ASSERT(m_pInfo);
	//I change colors less than 256 to grey as follows.
	if(GetBitsPerPixel()<=8)//if(m_bIsColTrue==false)
	{
		m_nPalEntries=GetPalEntries();
		RGBQUAD* pColTab=GetColorTablePtr();
		//Roll through the color table,and change each color to the relevant grey level.
		for(int ndx=0;ndx<m_nPalEntries;ndx++)
		{
			WORD wGrey=(WORD)pColTab[ndx].rgbRed*77+(WORD)pColTab[ndx].rgbGreen*150+(WORD)pColTab[ndx].rgbBlue*29;
			wGrey=wGrey/256;
			pColTab[ndx].rgbRed=(BYTE)wGrey;
			pColTab[ndx].rgbGreen=(BYTE)wGrey;
			pColTab[ndx].rgbBlue=(BYTE)wGrey;
			//pColTab[ndx].rgbReserved=0;
		}
	}
	else//I try to change True color to grey.
	{
		int npadw=WIDTH(GetWidth());
		int nNewPixels=GetHeight()*npadw;
		BYTE* pPixelTmp=new BYTE[nNewPixels];
		for(int ny=0;ny<GetHeight();ny++) for(int nx=0;nx<GetWidth();nx++)
		{
			BYTE* pbt=GetPixelByte(nx,ny);
			WORD wGrey=(WORD)pbt[2]*77+(WORD)pbt[1]*150+(WORD)pbt[0]*29;
			wGrey=wGrey/256; pPixelTmp[(GetHeight()-1-ny)*npadw+nx]=(BYTE)wGrey;
		}

		//change pixels here
		delete[] m_pPixels; m_pPixels=pPixelTmp; m_nPixels=nNewPixels;
		//Allocate space for a new bitmap info header,and copy the old info header.
		m_pInfo->bmiHeader.biBitCount=8;
		int nInfo=sizeof(BITMAPINFOHEADER)+1024;
		BITMAPINFO* pInfoTmp=(BITMAPINFO*)new BYTE[nInfo];
		memcpy(pInfoTmp,m_pInfo,sizeof(BITMAPINFOHEADER));
		delete[] (BYTE*)m_pInfo; m_pInfo=pInfoTmp;
		RGBQUAD* pColTab=GetColorTablePtr();
		m_nPalEntries=256;
		for(int ndx=0;ndx<m_nPalEntries;ndx++)
		{
			pColTab[ndx].rgbRed=(BYTE)ndx;
			pColTab[ndx].rgbGreen=(BYTE)ndx;
			pColTab[ndx].rgbBlue=(BYTE)ndx;
			pColTab[ndx].rgbReserved=0;
		}
		m_bIsColTrue=false;
	}
}

void CDIBitmap::SortGreyLevel()
{
	ASSERT(m_pPixels);
	if(GetBitsPerPixel()!=8) return;
	if(GetPalEntries()!=256) return;

	BYTE btGrey;
	int x,y;
	RGBQUAD* pRGB=GetColorTablePtr();
	DWORD dwPadWidth=WIDTH(GetWidth());

	for(y=0;y<GetHeight();y++) for(x=0;x<GetWidth();x++)
	{
		btGrey=m_pPixels[y*dwPadWidth+x];//btGrey=*GetPixelByte(x,y);
		btGrey=pRGB[(unsigned int)btGrey].rgbRed;
		m_pPixels[y*dwPadWidth+x]=btGrey;//SetPixelByte(x,y,&btGrey,1);
	}
	for(x=0;x<256;x++)
	{
		pRGB[x].rgbRed=(BYTE)x;
		pRGB[x].rgbGreen=(BYTE)x;
		pRGB[x].rgbBlue=(BYTE)x;
	}
}
