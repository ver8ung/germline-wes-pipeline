@echo off
REM Convenience wrapper so the pipeline runs by simply typing `run ...` from
REM cmd.exe OR PowerShell, without the .ps1 file-association opening Notepad and
REM without execution-policy friction. All arguments are forwarded to run.ps1.
REM
REM   run smoke -n          dry-run the smoke test
REM   run smoke             full smoke test
REM   run --cores 8         full pipeline
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0run.ps1" %*
