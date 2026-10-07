@echo off
chcp 65001 >nul
cd /d "%~dp0"
where python >nul 2>&1
if errorlevel 1 (
  echo 未找到 Python。请先安装 Python 3.9 或更高版本。
  exit /b 1
)
python -c "import sys; raise SystemExit(0 if sys.version_info >= (3, 9) else 1)"
if errorlevel 1 (
  echo 需要 Python 3.9 或更高版本。
  exit /b 1
)
python -m venv .venv
if errorlevel 1 exit /b 1
.venv\Scripts\python.exe -m pip install --upgrade pip
if errorlevel 1 exit /b 1
.venv\Scripts\python.exe -m pip install -r requirements.txt
if errorlevel 1 exit /b 1
if not exist .env copy /Y .env.example .env >nul
echo 环境已就绪: %cd%\.venv
echo 如需大模型，把密钥填进 .env 的 LLM_API_KEY。留空则只运行计算。
echo 编辑 input\request.txt 后运行 run.bat
