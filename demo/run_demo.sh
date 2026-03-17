#!/bin/bash
# 分子筛选工作流 Demo 运行脚本
# Molecular Screening Workflow Demo Script

echo "=========================================="
echo "  分子筛选工作流 - Demo 演示"
echo "=========================================="
echo ""

# 检查Python环境
if ! command -v python &> /dev/null; then
    echo "错误: 未找到Python，请先安装Python 3.8+"
    exit 1
fi

echo "Python版本:"
python --version
echo ""

# 检查依赖
echo "检查依赖包..."
python -c "import rdkit; print('✓ RDKit已安装')" || echo "✗ RDKit未安装"
python -c "import tensorflow; print('✓ TensorFlow已安装')" || echo "✗ TensorFlow未安装"
python -c "import pandas; print('✓ Pandas已安装')" || echo "✗ Pandas未安装"
echo ""

# 设置Demo模式
export DEMO_MODE=true

echo "开始运行Demo..."
echo "使用数据:"
echo "  - 训练数据: demo/data/sample_training.csv"
echo "  - 待筛选分子: demo/data/sample_molecules.csv"
echo ""

# 运行工作流
cd ..
python run_workflow.py --demo

echo ""
echo "=========================================="
echo "  Demo运行完成！"
echo "=========================================="
echo ""
echo "结果文件:"
ls -la results/*.csv 2>/dev/null || echo "  (结果文件将在results目录生成)"
echo ""
echo "查看详细日志: logs/workflow_*.log"
