@ECHO OFF
IF "%1"=="" GOTO Stop

echo Building Linkage Mapper Release number %1

@ECHO ON
xcopy /E/I/Y ..\demo ..\LinkageMapper"%1"\LM_demo
xcopy /Y ..\toolbox\doc\*.pdf ..\LinkageMapper"%1"\*.pdf
xcopy /E/I/Y ..\toolbox\scripts ..\LinkageMapper"%1"\toolbox\scripts
xcopy /E/I/Y ..\toolbox\styles ..\LinkageMapper"%1"\toolbox\styles
xcopy /Y ..\toolbox\*.tbx ..\LinkageMapper"%1"\toolbox\*.tbx

del ..\LinkageMapper"%1"\toolbox\scripts\cc*.*
del ..\LinkageMapper"%1"\toolbox\climate*.*
del ..\LinkageMapper"%1"\LM_demo\cc*.*
del ..\LinkageMapper"%1"\LM_demo\demoData\cc*.*

@ECHO OFF
if "%2"=="cc" GOTO CCCopy
if "%3"=="cc" GOTO CCCopy
if "%4"=="cc" GOTO CCCopy
GOTO CMTest

:CCCopy
@ECHO ON
xcopy /Y ..\toolbox\scripts\cc*.* ..\LinkageMapper"%1"\toolbox\scripts\cc*.*
xcopy /Y ..\toolbox\climate*.tbx ..\LinkageMapper"%1"\toolbox\climate*.tbx
xcopy /Y ..\demo\cc*.* ..\LinkageMapper"%1"\LM_demo\cc*.*
xcopy /Y ..\demo\demoData\cc*.* ..\LinkageMapper"%1"\LM_demo\demoData\cc*.*
xcopy /Y ..\toolbox\doc\climate*.pdf ..\LinkageMapper"%1"\climate*.pdf

:CMTest
REM @ECHO OFF
REM if "%2"=="cm" GOTO CMCopy
REM if "%3"=="cm" GOTO CMCopy
REM if "%4"=="cm" GOTO CMCopy
REM GOTO BMTest



:Finish
@ECHO OFF
echo releaseNum = "%1" > ..\LinkageMapper"%1"\toolbox\scripts\lm_version.py
echo.
echo Finished building Linkage Mapper release number %1
if "%2"=="" GOTO END
echo Options: %2 %3 %4
GOTO END

:Stop
@ECHO OFF
echo.
echo Error: must supply a version string
echo syntax: lm_build (version) (option to include climate tool)
echo example: lm_build 0.8.4 cc
echo.

:End