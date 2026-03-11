# 🧬 分子筛选工作流 (Molecular Screening Workflow)

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Active-success.svg)]()

> 基于RNN的分子生成与多阶段虚拟筛选工作流，用于药物发现中的先导化合物筛选

## 📖 目录

- [项目简介](#项目简介)
- [核心功能](#核心功能)
- [快速开始](#快速开始)
- [安装指南](#安装指南)
- [使用方法](#使用方法)
- [工作流架构](#工作流架构)
- [目录结构](#目录结构)
- [Demo演示](#demo演示)
- [结果说明](#结果说明)
- [常见问题](#常见问题)
- [贡献指南](#贡献指南)
- [许可证](#许可证)

## 🎯 项目简介

本项目是一个端到端的分子筛选工作流，整合了多种计算化学和机器学习技术，用于从大规模分子库中筛选潜在的药物候选分子。工作流采用模块化设计，每个阶段可独立运行，也可通过统一入口一键执行全流程。

### 核心特点

- ✅ **一键运行**: 单条命令完成全流程筛选
- ✅ **模块化设计**: 各阶段独立，易于定制和扩展
- ✅ **自动化流程**: 自动处理数据流转和依赖关系
- ✅ **详细日志**: 记录每个步骤的执行状态和耗时
- ✅ **异常处理**: 失败时保留中间结果，便于调试
- ✅ **可视化报告**: 自动生成评估图表和统计报告

## 🔬 核心功能

### 1. RNN分子生成
- 基于SELFIES表示的分子生成
- 使用LSTM网络学习分子结构
- 支持温度采样和Top-k采样
- 生成多样化、有效的分子

### 2. QSAR活性预测
- 基于机器学习的活性预测模型
- 支持多种分子描述符
- 自动特征选择和模型优化
- 输出IC50预测值

### 3. ADMET性质过滤
- 12项ADMET指标评估
- P-gp底物/抑制剂预测
- CYP酶相互作用预测
- 毒性评估（AMES、hERG、肝毒性等）

### 4. 分子对接
- 支持多种对接引擎（AutoDock Vina、GNINA、DiffDock）
- 批量对接处理
- 结合能评分和构象分析

### 5. SA Score评估
- 合成可及性评分
- 分子复杂度分析
- 筛选易于合成的候选分子

### 6. 结果可视化
- 化学空间分布图
- 活性预测散点图
- 性质分布直方图
- 多维度评估雷达图

## 🚀 快速开始

### 前置要求

- Python 3.8 或更高版本
- Conda（推荐）或 pip

### 一键安装

```bash
# 克隆仓库
git clone https://github.com/yourusername/molecular-screening-workflow.git
cd molecular-screening-workflow

# 创建虚拟环境
conda create -n mol-screen python=3.8
conda activate mol-screen

# 安装依赖
pip install -r requirements.txt

# 运行工作流
python run_workflow.py --demo
```

### 最小化运行

```bash
# 使用默认配置运行完整工作流
python run_workflow.py

# 查看帮助
python run_workflow.py --help
```

## 📦 安装指南

### 方法1: 使用Conda（推荐）

```bash
# 创建新环境
conda create -n mol-screen python=3.8
conda activate mol-screen

# 安装RDKit（必须先安装）
conda install -c conda-forge rdkit

# 安装其他依赖
pip install -r requirements.txt
```

### 方法2: 使用pip

```bash
# 创建虚拟环境
python -m venv venv
source venv/bin/activate  # Linux/Mac
# 或
venv\Scripts\activate  # Windows

# 安装依赖
pip install -r requirements.txt
```

### 验证安装

```bash
# 运行测试
python -m pytest tests/

# 运行Demo
python run_workflow.py --demo
```

## 💻 使用方法

### 基本用法

```bash
# 1. 使用默认配置运行
python run_workflow.py

# 2. 使用自定义配置文件
python run_workflow.py --config config/my_config.yaml

# 3. 运行演示模式
python run_workflow.py --demo

# 4. 只运行特定步骤
python run_workflow.py --steps rnn qsar admet
```

### 配置文件

编辑 `config/workflow_config.yaml` 来自定义工作流参数：

```yaml
# 数据路径
data:
  rnn_train_data: "data/your_molecules.csv"
  qsar_train_data: "data/your_activity_data.csv"

# 模型参数
rnn:
  num_molecules: 1000
  temperature: 0.7

qsar:
  ic50_threshold: 100
  top_n: 100
```

### 分步运行

如果需要单独运行某个模块：

```bash
# 只运行RNN生成
python src/rnn_workflow.py

# 只运行QSAR预测
python src/qsar_engine.py --input results/generated.csv

# 只运行ADMET过滤
python src/admet_engine.py --input results/qsar_output.csv
```

## 🏗️ 工作流架构

```
┌─────────────────────────────────────────────────────────────┐
│                    分子筛选工作流                              │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│  步骤1: RNN分子生成                                           │
│  - 输入: 训练分子库 (SMILES)                                  │
│  - 输出: 生成的分子 (SMILES)                                  │
│  - 技术: LSTM + SELFIES                                      │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│  步骤2: QSAR活性预测                                          │
│  - 输入: 生成的分子                                           │
│  - 输出: IC50预测值 + 筛选结果                                │
│  - 技术: Random Forest / Ensemble                            │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│  步骤3: ADMET性质过滤                                         │
│  - 输入: QSAR筛选后的分子                                     │
│  - 输出: ADMET合格的分子                                      │
│  - 指标: 12项ADMET性质                                        │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│  步骤4: 分子对接                                              │
│  - 输入: ADMET合格的分子 + 靶蛋白PDB                          │
│  - 输出: 对接构象 + 结合能                                    │
│  - 工具: AutoDock Vina / GNINA                               │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│  步骤5: SA Score评估                                          │
│  - 输入: 对接结果                                             │
│  - 输出: 合成可及性评分                                       │
│  - 目标: 筛选易于合成的分子                                   │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│  步骤6: 结果可视化                                            │
│  - 输出: 评估图表、统计报告                                   │
│  - 内容: 化学空间分布、活性分布、性质雷达图                    │
└─────────────────────────────────────────────────────────────┘
```

## 📁 目录结构

```
molecular-screening-workflow/
│
├── run_workflow.py          # 统一入口程序
├── requirements.txt          # 依赖包列表
├── .gitignore               # Git忽略文件
├── LICENSE                  # 许可证
├── README.md                # 项目说明文档
│
├── config/                  # 配置文件目录
│   └── workflow_config.yaml # 工作流配置
│
├── src/                     # 源代码目录
│   ├── rnn_workflow.py      # RNN分子生成
│   ├── qsar_engine.py       # QSAR预测引擎
│   ├── admet_engine.py      # ADMET过滤引擎
│   ├── submit_job_docking.py # 分子对接
│   ├── sa_scorer.py         # SA Score评估
│   └── visualizer.py        # 可视化模块
│
├── data/                    # 数据目录
│   ├── Anticancer_compounds@1994_from_paper.csv
│   └── ChEMBL_Cleaned_akt1-activities_zip.csv
│
├── models/                  # 模型文件目录
│   ├── selfies_generator_rnn.keras
│   └── ultimate_ensemble_qsar_model.pkl
│
├── results/                 # 结果输出目录
│   ├── QSAR/               # QSAR结果
│   ├── ADMET/              # ADMET结果
│   ├── docking/            # 对接结果
│   └── SA/                 # SA Score结果
│
├── logs/                    # 日志文件目录
│
├── demo/                    # 演示数据目录
│   ├── data/               # 最小数据集
│   ├── run_demo.sh         # Demo运行脚本
│   └── expected_result/    # 预期结果示例
│
└── docs/                    # 文档目录
    ├── workflow_architecture.md
    ├── module_description.md
    ├── result_format.md
    └── customization_guide.md
```

## 🎮 Demo演示

### 快速验证

```bash
# 进入Demo目录
cd demo

# 运行Demo脚本
bash run_demo.sh  # Linux/Mac
# 或
run_demo.bat      # Windows

# 查看结果
ls expected_result/
```

### Demo数据说明

Demo包含最小可运行数据集：
- `sample_training.csv`: 100个训练分子
- `sample_molecules.csv`: 50个待筛选分子
- `sample_receptor.pdb`: 示例靶蛋白

### 预期结果

运行Demo后，将在 `results/` 目录生成：
- `generated_smiles.csv`: 生成的分子
- `qsar_predictions.csv`: QSAR预测结果
- `admet_passed.csv`: ADMET过滤结果
- `docking_results.csv`: 对接结果
- `evaluation_summary.png`: 可视化报告

## 📊 结果说明

### 输出文件格式

#### 1. QSAR预测结果 (`top_100_generated_molecules_with_predictions.csv`)

| SMILES | pred_IC50 | pIC50 |
|--------|-----------|-------|
| CCO... | 45.23 | 7.34 |
| CCN... | 82.15 | 7.09 |

#### 2. ADMET过滤结果 (`admet_passed.csv`)

包含通过12项ADMET指标的分子及其性质

#### 3. 对接结果 (`docking_results.csv`)

| SMILES | binding_affinity | RMSD |
|--------|------------------|------|
| CCO... | -8.5 | 1.2 |

#### 4. SA Score结果 (`sa_filtered.csv`)

| SMILES | sa_score | rank |
|--------|----------|------|
| CCO... | 3.2 | 1 |

### 可视化报告

`evaluation_summary.png` 包含：
- QSAR预测 vs 实际活性散点图
- 化学空间分布图（MW vs LogP）
- IC50分布直方图
- SA Score箱线图

## ❓ 常见问题

### Q1: 安装RDKit失败怎么办？

**A**: RDKit推荐使用Conda安装：
```bash
conda install -c conda-forge rdkit
```

### Q2: 运行时内存不足？

**A**: 减少生成的分子数量或批大小：
```yaml
rnn:
  num_molecules: 100  # 减少到100
  batch_size: 16      # 减少批大小
```

### Q3: 如何使用自己的数据？

**A**: 准备CSV文件，包含以下列：
- RNN训练: `SMILES` 列
- QSAR训练: `SMILES`, `IC50` 列

然后修改配置文件中的数据路径。

### Q4: 如何添加新的筛选标准？

**A**: 编辑 `config/workflow_config.yaml` 中的参数，或修改对应模块的源代码。

### Q5: 工作流中断后如何继续？

**A**: 工作流会保留中间结果，可以单独运行后续步骤：
```bash
python run_workflow.py --steps qsar admet
```

## 🤝 贡献指南

欢迎贡献代码、报告问题或提出建议！

1. Fork 本仓库
2. 创建特性分支 (`git checkout -b feature/AmazingFeature`)
3. 提交更改 (`git commit -m 'Add some AmazingFeature'`)
4. 推送到分支 (`git push origin feature/AmazingFeature`)
5. 开启 Pull Request

## 📄 许可证

本项目采用 MIT 许可证 - 详见 [LICENSE](LICENSE) 文件

## 📧 联系方式

- 项目主页: https://github.com/yourusername/molecular-screening-workflow
- 问题反馈: https://github.com/yourusername/molecular-screening-workflow/issues
- 邮箱: your.email@example.com

## 🙏 致谢

感谢以下开源项目：
- [RDKit](https://www.rdkit.org/) - 化学信息学工具包
- [TensorFlow](https://www.tensorflow.org/) - 深度学习框架
- [AutoDock Vina](http://vina.scripps.edu/) - 分子对接工具
- [SELFIES](https://github.com/aspuru-guzik-group/selfies) - 分子表示方法

---

**⭐ 如果这个项目对您有帮助，请给一个Star！**
