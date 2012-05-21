@ECHO OFF
IF "%1"=="" GOTO Stop
echo Building Linkage Mapper Release number %1

@ECHO ON
xcopy /E/I/Y ..\demo ..\build\LM_demo
xcopy /Y ..\toolbox\doc\*.pdf ..\build\*.pdf
xcopy /E/I/Y ..\toolbox\scripts ..\build\toolbox\scripts
xcopy /E/I/Y ..\toolbox\styles ..\build\toolbox\styles
xcopy /Y ..\toolbox\*.tbx ..\build\toolbox\*.tbx
@ECHO OFF

echo releaseNum = "%1" > ..\build\toolbox\scripts\lm_version.py
echo Finished building Linkage Mapper Release number %1
GOTO END

:Stop
@ECHO OFF
echo Error: must supply a version string
echo syntax: build (version)
echo example: build 0.7.8
echo .

:End