@echo off
rem Runs every Stage 3 deck in its own folder.
rem   run_all.bat "C:\path\to\lsdyna.exe" 4
setlocal
set SOLVER=%~1
if "%SOLVER%"=="" set SOLVER=lsdyna
set NCPU=%~2
if "%NCPU%"=="" set NCPU=4

echo === stage3_lance_only ===
pushd "%~dp0stage3_lance_only"
"%SOLVER%" i=stage3_lance_only.k ncpu=%NCPU% memory=200m
popd

echo === stage3_extract_mu000 ===
pushd "%~dp0stage3_extract_mu000"
"%SOLVER%" i=stage3_extract_mu000.k ncpu=%NCPU% memory=200m
popd

echo === stage3_extract_mu020 ===
pushd "%~dp0stage3_extract_mu020"
"%SOLVER%" i=stage3_extract_mu020.k ncpu=%NCPU% memory=200m
popd

echo === stage3_extract_mu030 ===
pushd "%~dp0stage3_extract_mu030"
"%SOLVER%" i=stage3_extract_mu030.k ncpu=%NCPU% memory=200m
popd

echo === stage3_insert_mu020 ===
pushd "%~dp0stage3_insert_mu020"
"%SOLVER%" i=stage3_insert_mu020.k ncpu=%NCPU% memory=200m
popd

echo === stage3_extract_tpa_g010 ===
pushd "%~dp0stage3_extract_tpa_g010"
"%SOLVER%" i=stage3_extract_tpa_g010.k ncpu=%NCPU% memory=200m
popd

echo === stage3_extract_tpa_g085 ===
pushd "%~dp0stage3_extract_tpa_g085"
"%SOLVER%" i=stage3_extract_tpa_g085.k ncpu=%NCPU% memory=200m
popd

echo === stage3_extract_mu020_fine ===
pushd "%~dp0stage3_extract_mu020_fine"
"%SOLVER%" i=stage3_extract_mu020_fine.k ncpu=%NCPU% memory=200m
popd

echo === stage3_extract_mu020_slow ===
pushd "%~dp0stage3_extract_mu020_slow"
"%SOLVER%" i=stage3_extract_mu020_slow.k ncpu=%NCPU% memory=200m
popd

echo Done. Next:
echo   python scripts\extract_stage3.py "%~dp0"
echo   python scripts\postprocess_stage3.py
