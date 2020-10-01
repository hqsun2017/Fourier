# Microsoft Developer Studio Project File - Name="Fourier" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=Fourier - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "Fourier.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "Fourier.mak" CFG="Fourier - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Fourier - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "Fourier - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "Fourier - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GX /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x804 /d "NDEBUG"
# ADD RSC /l 0x804 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"..\Lib\Fourier.lib"

!ELSEIF  "$(CFG)" == "Fourier - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x804 /d "_DEBUG"
# ADD RSC /l 0x804 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"..\Lib\FourierD.lib"

!ENDIF 

# Begin Target

# Name "Fourier - Win32 Release"
# Name "Fourier - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\Color.cpp
# End Source File
# Begin Source File

SOURCE=.\DibMap.cpp
# End Source File
# Begin Source File

SOURCE=.\Filter.cpp
# End Source File
# Begin Source File

SOURCE=.\tFunc.cpp
# End Source File
# Begin Source File

SOURCE=.\x_Corr1D.cpp
# End Source File
# Begin Source File

SOURCE=.\xC2Base1D.cpp
# End Source File
# Begin Source File

SOURCE=.\xC2Base2D.cpp
# End Source File
# Begin Source File

SOURCE=.\xConv1D.cpp
# End Source File
# Begin Source File

SOURCE=.\xConv2D.cpp
# End Source File
# Begin Source File

SOURCE=.\xCorr1D.cpp
# End Source File
# Begin Source File

SOURCE=.\xCorr2D.cpp
# End Source File
# Begin Source File

SOURCE=.\xDCT_DST.cpp
# End Source File
# Begin Source File

SOURCE=.\xDeconv1D.cpp
# End Source File
# Begin Source File

SOURCE=.\xDeconv2D.cpp
# End Source File
# Begin Source File

SOURCE=.\xDiffusion2D.cpp
# End Source File
# Begin Source File

SOURCE=.\xFFT.cpp
# End Source File
# Begin Source File

SOURCE=.\xFFT2D.cpp
# End Source File
# Begin Source File

SOURCE=.\xFHT.cpp
# End Source File
# Begin Source File

SOURCE=.\xFHT2D.cpp
# End Source File
# Begin Source File

SOURCE=.\xHilbert.cpp
# End Source File
# Begin Source File

SOURCE=.\xHilbert2D.cpp
# End Source File
# Begin Source File

SOURCE=.\xPhase1D.cpp
# End Source File
# Begin Source File

SOURCE=.\xPhase2D.cpp
# End Source File
# Begin Source File

SOURCE=.\xPoisson2D.cpp
# End Source File
# Begin Source File

SOURCE=.\xRealFFT.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\Color.h
# End Source File
# Begin Source File

SOURCE=.\DibMap.h
# End Source File
# Begin Source File

SOURCE=.\Filter.h
# End Source File
# Begin Source File

SOURCE=.\tFunc.h
# End Source File
# Begin Source File

SOURCE=.\x_Corr1D.h
# End Source File
# Begin Source File

SOURCE=.\xC2Base1D.h
# End Source File
# Begin Source File

SOURCE=.\xC2Base2D.h
# End Source File
# Begin Source File

SOURCE=.\xConv1D.h
# End Source File
# Begin Source File

SOURCE=.\xConv2D.h
# End Source File
# Begin Source File

SOURCE=.\xCorr1D.h
# End Source File
# Begin Source File

SOURCE=.\xCorr2D.h
# End Source File
# Begin Source File

SOURCE=.\xDCT_DST.h
# End Source File
# Begin Source File

SOURCE=.\xDeconv1D.h
# End Source File
# Begin Source File

SOURCE=.\xDeconv2D.h
# End Source File
# Begin Source File

SOURCE=.\xDiffusion2D.h
# End Source File
# Begin Source File

SOURCE=.\xFFT.h
# End Source File
# Begin Source File

SOURCE=.\xFFT2D.h
# End Source File
# Begin Source File

SOURCE=.\xFHT.h
# End Source File
# Begin Source File

SOURCE=.\xFHT2D.h
# End Source File
# Begin Source File

SOURCE=.\xHilbert.h
# End Source File
# Begin Source File

SOURCE=.\xHilbert2D.h
# End Source File
# Begin Source File

SOURCE=.\xPhase1D.h
# End Source File
# Begin Source File

SOURCE=.\xPhase2D.h
# End Source File
# Begin Source File

SOURCE=.\xPoisson2D.h
# End Source File
# Begin Source File

SOURCE=.\xRealFFT.h
# End Source File
# End Group
# End Target
# End Project
