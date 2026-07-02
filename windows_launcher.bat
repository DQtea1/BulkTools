@echo off
setlocal enabledelayedexpansion

set "BROWSE_DIR=%USERPROFILE%"
set "OUT_DIR=%USERPROFILE%\shiny_out"

if not exist "%OUT_DIR%" mkdir "%OUT_DIR%"

rem --- Allocate ~3/4 of the logical CPU cores to the container ---
set /a CORES=%NUMBER_OF_PROCESSORS% * 3 / 4
if %CORES% LSS 1 set CORES=1
echo Using %CORES% of %NUMBER_OF_PROCESSORS% CPU core(s).

rem --- Mount every extra/external drive (D: onward) so it is browsable in the
rem     app as /mounts/<letter>. Drive C: is already covered by the home mount.
set "DRIVE_MOUNTS="
for %%D in (D E F G H I J K L M N O P Q R S T U V W X Y Z) do (
  if exist "%%D:\" (
    echo Mounting drive %%D: -^> /mounts/%%D
    set "DRIVE_MOUNTS=!DRIVE_MOUNTS! --mount type=bind,source=%%D:\,target=/mounts/%%D"
  )
)

start "" /b cmd /c "timeout /t 3 /nobreak >nul && start "" http://localhost:5288/"

docker run --rm -p 5288:5288 ^
  --cpus=%CORES% ^
  -e SHINY_PORT=5288 ^
  -e SHINY_ROOT_PATH=/browse ^
  -e SHINY_ROOT_NAME=home ^
  -e SHINY_N_CORES=%CORES% ^
  --mount type=bind,source="%BROWSE_DIR%",target=/browse ^
  --mount type=bind,source="%OUT_DIR%",target=/out ^
  !DRIVE_MOUNTS! ^
  qtea1/bulktools:latest

pause