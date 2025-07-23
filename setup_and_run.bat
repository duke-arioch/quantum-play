@echo off
echo Setting up Quantum Observer Dependence Test...
echo.

REM Check if Python is installed
python --version >nul 2>&1
if errorlevel 1 (
    echo Error: Python is not installed or not in PATH
    echo Please install Python 3.8 or higher from https://www.python.org/
    pause
    exit /b 1
)

REM Create virtual environment if it doesn't exist
if not exist "venv" (
    echo Creating virtual environment...
    python -m venv venv
    if errorlevel 1 (
        echo Error creating virtual environment
        pause
        exit /b 1
    )
)

REM Activate virtual environment and install dependencies
echo Activating virtual environment and installing dependencies...
call venv\Scripts\activate.bat
pip install -r requirements.txt

if errorlevel 1 (
    echo Error installing dependencies
    pause
    exit /b 1
)

echo.
echo Setup complete! Running the quantum observer test...
echo.
python qi2.py

echo.
echo Script completed. Press any key to exit...
pause >nul
