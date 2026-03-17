@echo off
REM 分子筛选工作流 Demo 运行脚本 (Windows)
REM Molecular Screening Workflow Demo Script

echo ==========================================
echo   分子筛选工作流 - Demo 演示
echo ==========================================
echo.

REM 检查Python环境
python --version >nul 2>&1
if errorlevel 1 (
    echo 错误: 未找到Python，请先安装Python 3.8+
    exit /b 1
)

echo Python版本:
python --version
echo.

REM 检查依赖
echo 检查依赖包...
python -c "import rdkit; print('✓ RDKit已安装')" 2>nul || echo ✗ RDKit未安装
python -c "import tensorflow; print('✓ TensorFlow已安装')" 2>nul || echo ✗ TensorFlow未安装
python -c "import pandas; print('✓ Pandas已安装')" 2>nul || echo ✗ Pandas未安装
echo.

REM 设置Demo模式
set DEMO_MODE=true

echo 开始运行Demo...
echo 使用数据:
echo   - 训练数据: demo/data/sample_training.csv
echo   - 待筛选分子: demo/data/sample_molecules.csv
echo.

REM 运行工作流
cd ..
python run_workflow.py --demo

echo.
echo ==========================================
echo   Demo运行完成！
echo ==========================================
echo.
echo 结果文件:
dir results\*.csv 2>nul || echo   (结果文件将在results目录生成)
echo.
echo 查看详细日志: logs\workflow_*.log

pause
