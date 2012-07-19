@ECHO OFF
IF "%1"=="" GOTO Stop

echo Building Linkage Mapper Release number %1

@ECHO ON
xcopy /E/I/Y ..\demo ..\build\LM_demo
xcopy /Y ..\toolbox\doc\*.pdf ..\build\*.pdf
xcopy /E/I/Y ..\toolbox\scripts ..\build\toolbox\scripts
xcopy /E/I/Y ..\toolbox\styles ..\build\toolbox\styles
xcopy /Y ..\toolbox\*.tbx ..\build\toolbox\*.tbx

del ..\build\toolbox\scripts\cc*.*
del ..\build\toolbox\scripts\s7*.*
del ..\build\toolbox\scripts\s6*.*
del ..\build\toolbox\scripts\bar*.*
del ..\build\toolbox\climate*.*
del ..\build\LM_demo\cc*.*
del ..\build\LM_demo\demoData\cc*.*
del ..\build\toolbox\*Plus.tbx
del ..\build\toolbox\*Extras.tbx

@ECHO OFF
if "%2"=="cc" GOTO CCCopy
if "%3"=="cc" GOTO CCCopy
if "%4"=="cc" GOTO CCCopy
GOTO CMTest

:CCCopy
@ECHO ON
xcopy /Y ..\toolbox\scripts\cc*.* ..\build\toolbox\scripts\cc*.*
xcopy /Y ..\toolbox\climate*.tbx ..\build\toolbox\climate*.tbx
xcopy /Y ..\demo\cc*.* ..\build\LM_demo\cc*.*
xcopy /Y ..\demo\demoData\cc*.* ..\build\LM_demo\demoData\cc*.*

:CMTest
@ECHO OFF
if "%2"=="cm" GOTO CMCopy
if "%3"=="cm" GOTO CMCopy
if "%4"=="cm" GOTO CMCopy
GOTO BMTest

:CMCopy
@ECHO ON
xcopy /Y ..\toolbox\scripts\s7_c*.* ..\build\toolbox\scripts\s7_c*.*
xcopy /Y ..\toolbox\*Extras.tbx ..\build\toolbox\*Extras.tbx

:BMTest
@ECHO OFF
if "%2"=="bm" GOTO BMCopy
if "%3"=="bm" GOTO BMCopy
if "%4"=="bm" GOTO BMCopy
GOTO Finish

:BMCopy
@ECHO ON
xcopy /Y ..\toolbox\scripts\bar*.* ..\build\toolbox\scripts\bar*.*
xcopy /Y ..\toolbox\scripts\s6_b*.* ..\build\toolbox\scripts\s6_b*.*
xcopy /Y ..\toolbox\*Extras.tbx ..\build\toolbox\*Extras.tbx

:Finish
@ECHO OFF
echo releaseNum = "%1" > ..\build\toolbox\scripts\lm_version.py
echo.
echo Finished building Linkage Mapper release number %1
if "%2"=="" GOTO END
echo Options: %2 %3 %4
GOTO END

:Stop
@ECHO OFF
echo Error: must supply a version string
echo syntax: lm_build (version) cm bm cc
echo example: lm_build 0.8.4 cm cc
echo.

:End