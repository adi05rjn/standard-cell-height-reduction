@echo off
rem Builds a Manager project:   build_project.bat <variant> [nobuild]
rem   variant: 6.5 | 5.5 | 4.5 | handshake | bpr
rem 1) writes layout, layer map, model file and process JSON  2) converts the JSON to a .vproc
rem 3) links the material file  4) builds the model with ModelBuilder, as the Manager would
setlocal
call "C:\Coventor\SEMulator3D13\bin\SEMulator3D.bat" >nul 2>&1
rem COVLD_LICENSE_FILE (port@host of your SEMulator3D license server) must be set in your environment before running.
cd /d "%~dp0"
set T=%1
for /f %%i in ('python make_project.py %T% --name') do set N=%%i
if exist "%~dp0%N%\%N%-Output" rmdir /s /q "%~dp0%N%\%N%-Output"
if exist "%~dp0%N%\%N%-Input" rmdir /s /q "%~dp0%N%\%N%-Input"
python make_project.py %T% || exit /b 1
cd /d "%~dp0%N%"
if exist %N%.vproc del %N%.vproc
(echo %N%.json&echo %N%.vproc&echo y) | python "C:\Coventor\SEMulator3D13\apps\contrib\py\JSONTranslator.pyc" --noComments
rem the translator names its output "<name>.vproc .vproc"; rename it
if exist "%N%.vproc .vproc" ren "%N%.vproc .vproc" %N%.vproc
cd /d "%~dp0"
python link_materials.py %N% %N% || exit /b 1
if "%2"=="nobuild" exit /b 0
cd /d "%~dp0%N%"
ModelBuilder.exe -m %N%.zam
echo exit=%errorlevel%
