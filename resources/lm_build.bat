@ECHO OFF
IF "%1"=="" GOTO Stop

echo Building Linkage Mapper Release number %1

SET "rel_text=%1"
SETLOCAL EnableDelayedExpansion
SET rel_text=!rel_text:.=_!

@ECHO ON
xcopy /E/I/Y ..\demo ..\LinkageMapper"%rel_text%"\LM_demo
xcopy /E/I/Y ..\demo ..\LinkageMapper"%rel_text%"\LM_demo
xcopy /Y ..\toolbox\doc\*.pdf ..\LinkageMapper"%rel_text%"\*.pdf
xcopy /E/I/Y ..\toolbox\scripts ..\LinkageMapper"%rel_text%"\toolbox\scripts
xcopy /E/I/Y ..\toolbox\styles ..\LinkageMapper"%rel_text%"\toolbox\styles
xcopy /Y ..\toolbox\*.tbx ..\LinkageMapper"%rel_text%"\toolbox\*.tbx

del ..\LinkageMapper"%rel_text%"\toolbox\scripts\cc*.*
del ..\LinkageMapper"%rel_text%"\toolbox\climate*.*
del ..\LinkageMapper"%rel_text%"\LM_demo\cc*.*
del ..\LinkageMapper"%rel_text%"\LM_demo\demoData\cc*.*
del ..\LinkageMapper"%rel_text%"\toolbox\scripts\iterate*.*

REM @ECHO OFF
REM if "%2"=="cc" GOTO CCCopy
REM if "%3"=="cc" GOTO CCCopy
REM if "%4"=="cc" GOTO CCCopy
REM GOTO CMTest

:CCCopy
@ECHO ON
xcopy /Y ..\toolbox\scripts\cc*.* ..\LinkageMapper"%rel_text%"\toolbox\scripts\cc*.*
xcopy /Y ..\toolbox\climate*.tbx ..\LinkageMapper"%rel_text%"\toolbox\climate*.tbx
xcopy /Y ..\demo\cc*.* ..\LinkageMapper"%rel_text%"\LM_demo\cc*.*
xcopy /Y ..\demo\demoData\cc*.* ..\LinkageMapper"%rel_text%"\LM_demo\demoData\cc*.*
xcopy /Y ..\toolbox\doc\climate*.pdf ..\LinkageMapper"%rel_text%"\climate*.pdf

:CMTest
REM @ECHO OFF
REM if "%2"=="cm" GOTO CMCopy
REM if "%3"=="cm" GOTO CMCopy
REM if "%4"=="cm" GOTO CMCopy
REM GOTO BMTest



:Finish
@ECHO OFF
echo releaseNum = "%1" > ..\LinkageMapper"%rel_text%"\toolbox\scripts\lm_version.py
echo.
echo Finished building Linkage Mapper release number %1
if "%2"=="" GOTO END
echo Options: %2 %3 %4
GOTO END
ENDLOCAL

:Stop
@ECHO OFF
echo.
echo Error: must supply a version string
echo syntax: lm_build (version)
echo example: lm_build 0.8.4
echo.

:End