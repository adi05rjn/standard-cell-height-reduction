@echo off
call "C:\Coventor\SEMulator3D13\bin\SEMulator3D.bat" >nul 2>&1
rem COVLD_LICENSE_FILE (port@host of your SEMulator3D license server) must be set in your environment before running.
cd /d "%~dp0"
python build_stage2.py %1 %2
echo exit=%errorlevel%
