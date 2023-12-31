@echo off

echo Converts all contents of folder 'input' to .nii. If multiple series are within a subfolder, they recveive a different suffix (_e1, _e2, ...). If a subfolder contains a time series, the time series will be saved to the same .nii output file. 

del /q /f output\*
cls


for /d %%i in ("input\*") do (
    echo --------------------------
    echo %%~nxi

    call MRIcroGL\Resources\dcm2niix.exe -f "%%~nxi" -p y -z n -b n -o "output" "input\%%~nxi"
)


pause
exit
