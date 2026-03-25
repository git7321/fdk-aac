# Microsoft Developer Studio Project File - Name="fdkaac" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103
# TARGTYPE "Win32 (x86) Dynamic-Link Library" 0x0102
# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=fdkaac - Win32 Release
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "FDKAAC.MAK".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "FDKAAC.MAK" CFG="fdkaac - Win32 Release"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "fdkaac - Win32 DLL Release" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE "fdkaac - Win32 LIB Release" (based on "Win32 (x86) Static Library")
!MESSAGE "fdkaac - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "..\"
# PROP BASE Intermediate_Dir "..\"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\"
# PROP Intermediate_Dir "..\"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
LIB32=link.exe
CPP=cl.exe
# ADD BASE CPP /nologo /MD /W3 /Gy /I "." /I ".\libFDK\include" /I ".\libAACenc\include" /I ".\libSYS\include" /I ".\libAACdec\include" /I ".\libArithCoding\include" /I ".\libMpegTPDec\include" /I ".\libMpegTPEnc\include" /I ".\libSBRdec\include" /I ".\libSBRenc\include" /I ".\libPCMutils\include" /I ".\libDRCdec\include" /I ".\libSACenc\include" /I ".\libSACdec\include" /D "NDEBUG" /D "_WINDOWS" /D "WIN32" /D "_MBCS" /D "FDK_DLL" /D "FDK_EXPORT" /FD /Oxs /GF /c
# ADD CPP /nologo /MD /W3 /Gy /I "." /I ".\libFDK\include" /I ".\libAACenc\include" /I ".\libSYS\include" /I ".\libAACdec\include" /I ".\libArithCoding\include" /I ".\libMpegTPDec\include" /I ".\libMpegTPEnc\include" /I ".\libSBRdec\include" /I ".\libSBRenc\include" /I ".\libPCMutils\include" /I ".\libDRCdec\include" /I ".\libSACenc\include" /I ".\libSACdec\include" /D "NDEBUG" /D "_WINDOWS" /D "WIN32" /D "_MBCS" /D "FDK_DLL" /D "FDK_EXPORT" /FD /Oxs /GF /c
MTL=midl.exe
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32 /win32
RSC=rc.exe
# ADD BASE RSC /l 0x409 /d "FDKAAC_DLL" /d "NDEBUG" /c 1252
# ADD RSC /l 0x409 /d "FDKAAC_DLL" /d "NDEBUG" /c 1252
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib /nologo /base:"0x79680000" /version:5.2 /subsystem:windows /dll /debug /machine:I386 /osversion:5.2 /opt:ref,icf,nowin98 /release
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib /nologo /base:"0x79680000" /version:5.2 /subsystem:windows /dll /debug /machine:I386 /osversion:5.2 /opt:ref,icf,nowin98 /release

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "..\"
# PROP BASE Intermediate_Dir "..\"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\"
# PROP Intermediate_Dir "..\"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
CPP=cl.exe
# ADD BASE CPP /nologo /MD /W3 /Gy /I "." /I ".\libFDK\include" /I ".\libAACenc\include" /I ".\libSYS\include" /I ".\libAACdec\include" /I ".\libArithCoding\include" /I ".\libMpegTPDec\include" /I ".\libMpegTPEnc\include" /I ".\libSBRdec\include" /I ".\libSBRenc\include" /I ".\libPCMutils\include" /I ".\libDRCdec\include" /I ".\libSACenc\include" /I ".\libSACdec\include" /D "NDEBUG" /D "_WINDOWS" /D "WIN32" /D "_MBCS" /D "_LIB" /FD /Oxs /GF /GL /c
# ADD CPP /nologo /MD /W3 /Gy /I "." /I ".\libFDK\include" /I ".\libAACenc\include" /I ".\libSYS\include" /I ".\libAACdec\include" /I ".\libArithCoding\include" /I ".\libMpegTPDec\include" /I ".\libMpegTPEnc\include" /I ".\libSBRdec\include" /I ".\libSBRenc\include" /I ".\libPCMutils\include" /I ".\libDRCdec\include" /I ".\libSACenc\include" /I ".\libSACdec\include" /D "NDEBUG" /D "_WINDOWS" /D "WIN32" /D "_MBCS" /D "_LIB" /FD /Oxs /GF /GL /c
RSC=rc.exe
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo /ltcg
# ADD LIB32 /nologo /ltcg

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "..\"
# PROP BASE Intermediate_Dir "..\"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\"
# PROP Intermediate_Dir "..\"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
MTL=midl.exe
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32
CPP=cl.exe
# ADD BASE CPP /nologo /MD /W3 /Gy /I "." /I ".\libSYS\include" /I ".\libAACenc\include" /D "NDEBUG" /D "_CONSOLE" /D "WIN32" /D "_MBCS" /FD /Oxs /GF /GL /c
# ADD CPP /nologo /MD /W3 /Gy /I "." /I ".\libSYS\include" /I ".\libAACenc\include" /D "NDEBUG" /D "_CONSOLE" /D "WIN32" /D "_MBCS" /FD /Oxs /GF /GL /c
RSC=rc.exe
# ADD BASE RSC /l 0x409 /d "NDEBUG" /c 1252
# ADD RSC /l 0x409 /d "NDEBUG" /c 1252
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib ..\fdkaac.lib uuid.lib /nologo /base:"0x796A0000" /version:5.2 /subsystem:console /debug /machine:I386 /osversion:5.2 /opt:ref,icf,nowin98 /release /ltcg
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib ..\fdkaac.lib uuid.lib /nologo /base:"0x796A0000" /version:5.2 /subsystem:console /debug /machine:I386 /osversion:5.2 /opt:ref,icf,nowin98 /release /ltcg

!ENDIF 

# Begin Target

# Name "fdkaac - Win32 DLL Release"
# Name "fdkaac - Win32 LIB Release"
# Name "fdkaac - Win32 Release"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Group "fdkaac"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\aacenc.c

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\aacenc.h
# End Source File
# Begin Source File

SOURCE=.\caf_reader.c

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\catypes.h
# End Source File
# Begin Source File

SOURCE=.\compat.h
# End Source File
# Begin Source File

SOURCE=.\compat_win32.c

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\extrapolater.c

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\getopt.c

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\getopt.h
# End Source File
# Begin Source File

SOURCE=.\limiter.c

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\lpc.c

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\lpc.h
# End Source File
# Begin Source File

SOURCE=.\lpcm.h
# End Source File
# Begin Source File

SOURCE=.\m4af.c

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\m4af.h
# End Source File
# Begin Source File

SOURCE=.\m4af_endian.h
# End Source File
# Begin Source File

SOURCE=.\main.c

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\metadata.c

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\metadata.h
# End Source File
# Begin Source File

SOURCE=.\parson.c

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\parson.h
# End Source File
# Begin Source File

SOURCE=.\pcm_float_converter.c

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\pcm_native_converter.c

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\pcm_reader.h
# End Source File
# Begin Source File

SOURCE=.\pcm_readhelper.c

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\pcm_sint16_converter.c

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\progress.c

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\progress.h
# End Source File
# Begin Source File

SOURCE=.\wav_reader.c

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

!ENDIF 

# End Source File
# End Group
# Begin Source File

SOURCE=.\libAACdec\src\aac_ram.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\aac_ram.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\aac_rom.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\aac_rom.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\aacdec_drc.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\aacdec_drc.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\aacdec_drc_types.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\aacdec_hcr.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\aacdec_hcr.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\aacdec_hcr_bit.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\aacdec_hcr_bit.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\aacdec_hcr_types.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\aacdec_hcrs.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\aacdec_hcrs.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\aacdec_pns.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\aacdec_pns.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\aacdec_tns.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\aacdec_tns.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\aacdecoder.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\aacdecoder.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\aacdecoder_lib.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\include\aacdecoder_lib.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\aacenc.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\aacenc.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\aacenc_lib.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\include\aacenc_lib.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\aacenc_pns.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\aacenc_pns.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\aacEnc_ram.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\aacEnc_ram.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\aacEnc_rom.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\aacEnc_rom.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\aacenc_tns.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\aacenc_tns.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\include\abs.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\include\x86\abs_x86.h
# End Source File
# Begin Source File

SOURCE=.\libArithCoding\src\ac_arith_coder.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libArithCoding\include\ac_arith_coder.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\adj_thr.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\adj_thr.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\adj_thr_data.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\src\autocorr2nd.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libFDK\include\autocorr2nd.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\band_nrg.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\band_nrg.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\bandwidth.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\bandwidth.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\bit_cnt.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\bit_cnt.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\bit_sbr.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\bit_sbr.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\bitenc.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\bitenc.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\block.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\block.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\block_switch.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\block_switch.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\channel.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\channel.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\channel_map.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\channel_map.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\channelinfo.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\channelinfo.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\chaosmeasure.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\chaosmeasure.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\include\clz.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\include\x86\clz_x86.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\cmondata.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\code_env.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\code_env.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\include\common_fix.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\conceal.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\conceal.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\conceal_types.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\include\cplx_mul.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\src\dct.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libFDK\include\dct.h
# End Source File
# Begin Source File

SOURCE=.\libDRCdec\src\drcDec_gainDecoder.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libDRCdec\src\drcDec_gainDecoder.h
# End Source File
# Begin Source File

SOURCE=.\libDRCdec\src\drcDec_reader.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libDRCdec\src\drcDec_reader.h
# End Source File
# Begin Source File

SOURCE=.\libDRCdec\src\drcDec_rom.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libDRCdec\src\drcDec_rom.h
# End Source File
# Begin Source File

SOURCE=.\libDRCdec\src\drcDec_selectionProcess.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libDRCdec\src\drcDec_selectionProcess.h
# End Source File
# Begin Source File

SOURCE=.\libDRCdec\src\drcDec_tools.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libDRCdec\src\drcDec_tools.h
# End Source File
# Begin Source File

SOURCE=.\libDRCdec\src\drcDec_types.h
# End Source File
# Begin Source File

SOURCE=.\libDRCdec\src\drcDecoder.h
# End Source File
# Begin Source File

SOURCE=.\libDRCdec\src\drcGainDec_init.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libDRCdec\src\drcGainDec_init.h
# End Source File
# Begin Source File

SOURCE=.\libDRCdec\src\drcGainDec_preprocess.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libDRCdec\src\drcGainDec_preprocess.h
# End Source File
# Begin Source File

SOURCE=.\libDRCdec\src\drcGainDec_process.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libDRCdec\src\drcGainDec_process.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\dyn_bits.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\dyn_bits.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\env_bit.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\env_bit.h
# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\env_calc.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\env_calc.h
# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\env_dec.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\env_dec.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\env_est.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\env_est.h
# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\env_extr.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\env_extr.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\include\FDK_archdef.h
# End Source File
# Begin Source File

SOURCE=.\libSYS\include\FDK_audio.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\src\FDK_bitbuffer.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libFDK\include\FDK_bitbuffer.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\include\FDK_bitstream.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\src\FDK_core.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libFDK\include\FDK_core.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\src\FDK_crc.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libFDK\include\FDK_crc.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\src\FDK_decorrelate.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libFDK\include\FDK_decorrelate.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\FDK_delay.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\FDK_delay.h
# End Source File
# Begin Source File

SOURCE=.\libDRCdec\src\FDK_drcDecLib.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libDRCdec\include\FDK_drcDecLib.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\src\FDK_hybrid.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libFDK\include\FDK_hybrid.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\src\FDK_lpc.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libFDK\include\FDK_lpc.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\src\FDK_matrixCalloc.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libFDK\include\FDK_matrixCalloc.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\src\FDK_qmf_domain.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libFDK\include\FDK_qmf_domain.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\src\FDK_tools_rom.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libFDK\include\FDK_tools_rom.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\src\FDK_trigFcts.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libFDK\include\FDK_trigFcts.h
# End Source File
# Begin Source File

SOURCE=".\fdkaac.def"

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libFDK\src\fft.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libFDK\include\fft.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\src\fft_rad2.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libFDK\include\fft_rad2.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\include\fixmadd.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\include\fixminmax.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\include\fixmul.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\include\x86\fixmul_x86.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\src\fixpoint_math.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libFDK\include\fixpoint_math.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\include\x86\fixpoint_math_x86.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\fram_gen.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\fram_gen.h
# End Source File
# Begin Source File

SOURCE=.\libSYS\src\genericStds.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSYS\include\genericStds.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\grp_data.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\grp_data.h
# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\hbe.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\hbe.h
# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\HFgen_preFlat.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\HFgen_preFlat.h
# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\huff_dec.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\huff_dec.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\src\huff_nodes.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libFDK\include\huff_nodes.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\intensity.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\intensity.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\interface.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\invf_est.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\invf_est.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\ldfiltbank.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\ldfiltbank.h
# End Source File
# Begin Source File

SOURCE=.\libPCMutils\src\limiter.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libPCMutils\include\limiter.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\line_pe.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\line_pe.h
# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\lpp_tran.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\lpp_tran.h
# End Source File
# Begin Source File

SOURCE=.\libSYS\include\machine_type.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\src\mdct.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libFDK\include\mdct.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\metadata_compressor.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\metadata_compressor.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\metadata_main.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\metadata_main.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\mh_det.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\mh_det.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\mps_main.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\mps_main.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\ms_stereo.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\ms_stereo.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\nf_est.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\nf_est.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\src\nlc_dec.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libFDK\include\nlc_dec.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\noisedet.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\noisedet.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\overlapadd.h
# End Source File
# Begin Source File

SOURCE=.\libPCMutils\src\pcm_utils.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libPCMutils\include\pcm_utils.h
# End Source File
# Begin Source File

SOURCE=.\libPCMutils\src\pcmdmx_lib.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libPCMutils\include\pcmdmx_lib.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\pns_func.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\pnsparam.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\pnsparam.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\pre_echo_control.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\pre_echo_control.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\ps_bitenc.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\ps_bitenc.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\ps_const.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\ps_encode.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\ps_encode.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\ps_main.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\ps_main.h
# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\psbitdec.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\psbitdec.h
# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\psdec.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\psdec.h
# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\psdec_drm.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\psdec_drm.h
# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\psdecrom_drm.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\psy_configuration.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\psy_configuration.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\psy_const.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\psy_data.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\psy_main.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\psy_main.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\pulsedata.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\pulsedata.h
# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\pvc_dec.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\pvc_dec.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\qc_data.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\qc_main.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\qc_main.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\src\qmf.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libFDK\include\qmf.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\include\qmf_pcm.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\quantize.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\quantize.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\resampler.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\resampler.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\rvlc.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\rvlc.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\rvlc_info.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\rvlcbit.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\rvlcbit.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\rvlcconceal.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\rvlcconceal.h
# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_bitdec.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_bitdec.h
# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_calcM1andM2.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_calcM1andM2.h
# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_dec.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_dec.h
# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_dec_conceal.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_dec_conceal.h
# End Source File
# Begin Source File

SOURCE=.\libSACdec\include\sac_dec_errorcodes.h
# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_dec_interface.h
# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_dec_lib.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACdec\include\sac_dec_lib.h
# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_dec_ssc_struct.h
# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_process.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_process.h
# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_qmf.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_qmf.h
# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_reshapeBBEnv.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_reshapeBBEnv.h
# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_rom.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_rom.h
# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_smoothing.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_smoothing.h
# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_stp.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_stp.h
# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_tsd.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACdec\src\sac_tsd.h
# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_bitstream.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_bitstream.h
# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_const.h
# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_delay.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_delay.h
# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_dmx_tdom_enh.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_dmx_tdom_enh.h
# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_filter.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_filter.h
# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_framewindowing.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_framewindowing.h
# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_huff_tab.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_huff_tab.h
# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_lib.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACenc\include\sacenc_lib.h
# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_nlc_enc.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_nlc_enc.h
# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_onsetdetect.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_onsetdetect.h
# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_paramextract.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_paramextract.h
# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_staticgain.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_staticgain.h
# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_tree.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_tree.h
# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_vectorfunctions.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSACenc\src\sacenc_vectorfunctions.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\sbr.h
# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\sbr_deb.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\sbr_deb.h
# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\sbr_dec.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\sbr_dec.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\sbr_def.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\sbr_encoder.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRenc\include\sbr_encoder.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\sbr_misc.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\sbr_misc.h
# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\sbr_ram.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\sbr_ram.h
# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\sbr_rom.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\sbr_rom.h
# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\sbrdec_drc.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\sbrdec_drc.h
# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\sbrdec_freq_sca.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\sbrdec_freq_sca.h
# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\sbrdecoder.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRdec\include\sbrdecoder.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\sbrenc_freq_sca.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\sbrenc_freq_sca.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\sbrenc_ram.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\sbrenc_ram.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\sbrenc_rom.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\sbrenc_rom.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\src\scale.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libFDK\include\scale.h
# End Source File
# Begin Source File

SOURCE=.\libFDK\include\scramble.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\sf_estim.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\sf_estim.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\spreading.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\spreading.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\stereo.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\stereo.h
# End Source File
# Begin Source File

SOURCE=.\libSYS\src\syslib_channelMapDescr.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSYS\include\syslib_channelMapDescr.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\tns_func.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\ton_corr.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\ton_corr.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\tonality.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\tonality.h
# End Source File
# Begin Source File

SOURCE=.\libMpegTPDec\include\tp_data.h
# End Source File
# Begin Source File

SOURCE=.\libMpegTPEnc\include\tp_data.h
# End Source File
# Begin Source File

SOURCE=.\libMpegTPDec\src\tp_version.h
# End Source File
# Begin Source File

SOURCE=.\libMpegTPEnc\src\tp_version.h
# End Source File
# Begin Source File

SOURCE=.\libMpegTPDec\src\tpdec_adif.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libMpegTPDec\src\tpdec_adif.h
# End Source File
# Begin Source File

SOURCE=.\libMpegTPDec\src\tpdec_adts.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libMpegTPDec\src\tpdec_adts.h
# End Source File
# Begin Source File

SOURCE=.\libMpegTPDec\src\tpdec_asc.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libMpegTPDec\src\tpdec_drm.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libMpegTPDec\src\tpdec_drm.h
# End Source File
# Begin Source File

SOURCE=.\libMpegTPDec\src\tpdec_latm.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libMpegTPDec\src\tpdec_latm.h
# End Source File
# Begin Source File

SOURCE=.\libMpegTPDec\src\tpdec_lib.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libMpegTPDec\include\tpdec_lib.h
# End Source File
# Begin Source File

SOURCE=.\libMpegTPEnc\src\tpenc_adif.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libMpegTPEnc\src\tpenc_adif.h
# End Source File
# Begin Source File

SOURCE=.\libMpegTPEnc\src\tpenc_adts.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libMpegTPEnc\src\tpenc_adts.h
# End Source File
# Begin Source File

SOURCE=.\libMpegTPEnc\src\tpenc_asc.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libMpegTPEnc\src\tpenc_asc.h
# End Source File
# Begin Source File

SOURCE=.\libMpegTPEnc\src\tpenc_latm.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libMpegTPEnc\src\tpenc_latm.h
# End Source File
# Begin Source File

SOURCE=.\libMpegTPEnc\src\tpenc_lib.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libMpegTPEnc\include\tpenc_lib.h
# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\tran_det.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libSBRenc\src\tran_det.h
# End Source File
# Begin Source File

SOURCE=.\libSBRdec\src\transcendent.h
# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\transform.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACenc\src\transform.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\usacdec_ace_d4t64.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\usacdec_ace_d4t64.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\usacdec_ace_ltp.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\usacdec_ace_ltp.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\usacdec_acelp.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\usacdec_acelp.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\usacdec_const.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\usacdec_fac.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\usacdec_fac.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\usacdec_lpc.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\usacdec_lpc.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\usacdec_lpd.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\usacdec_lpd.h
# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\usacdec_rom.cpp

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\libAACdec\src\usacdec_rom.h
# End Source File
# Begin Source File

SOURCE=.\libPCMutils\src\version.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# Begin Source File

SOURCE=.\fdkaac.rc

!IF  "$(CFG)" == "fdkaac - Win32 DLL Release"

!ELSEIF  "$(CFG)" == "fdkaac - Win32 LIB Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "fdkaac - Win32 Release"

!ENDIF 

# End Source File
# End Group
# End Target
# End Project
