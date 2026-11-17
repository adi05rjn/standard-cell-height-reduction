@echo off
rem Whole-model verification:  run_verify.bat <variant> [resolution]     variant: 6.5 | 5.5 | 4.5 | handshake | bpr
call "C:\Coventor\SEMulator3D13\bin\SEMulator3D.bat" >nul 2>&1
rem COVLD_LICENSE_FILE (port@host of your SEMulator3D license server) must be set in your environment before running.
cd /d "%~dp0"
python verify_model.py %1 %2
echo exit=%errorlevel%
