# 模块详细说明

本文档详细描述工作流中各个模块的功能、输入、输出和使用方法。

---

## 目录

1. [RNN分子生成模块](#rnn分子生成模块)
2. [QSAR活性预测模块](#qsar活性预测模块)
3. [ADMET过滤模块](#admet过滤模块)
4. [分子对接模块](#分子对接模块)
5. [SA Score评估模块](#sa-score评估模块)
6. [可视化模块](#可视化模块)
7. [日志模块](#日志模块)
8. [工具模块](#工具模块)

---

## RNN分子生成模块

### 文件位置
`src/rnn_workflow.py`

### 功能描述
使用循环神经网络（RNN）生成新的分子结构。基于SELFIES表示法，训练LSTM网络学习分子结构规律，然后生成多样化的新分子。

### 核心类/函数

#### `RNNMoleculeGenerator`

**初始化参数**:
```python
def __init__(self, config):
    self.num_molecules = config.get('num_molecules', 500)
    self.temperature = config.get('temperature', 0.7)
    self.top_k = config.get('top_k', 50)
    self.epochs = config.get('epochs', 100)
    self.batch_size = config.get('batch_size', 32)
```

**主要方法**:

1. `train(data_path)`
   - 功能: 训练RNN模型
   - 输入: 训练数据CSV文件路径
   - 输出: 训练好的模型文件

2. `generate(num_molecules)`
   - 功能: 生成新分子
   - 输入: 生成分子数量
   - 输出: SMILES列表

3. `save_model(path)`
   - 功能: 保存模型
   - 输入: 保存路径

4. `load_model(path)`
   - 功能: 加载模型
   - 输入: 模型路径

### 输入格式

**训练数据** (`data/Anticancer_compounds@1994_from_paper.csv`):
```csv
SMILES
CCOc1ccc(CC(=O)O)cc1
CC(C)Cc1ccc(C(C)C(=O)O)cc1
...
```

### 输出格式

**生成结果** (`results/generated_smiles.csv`):
```csv
SMILES
CCOc1ccc(CC(=O)O)cc1
CC(C)Cc1ccc(C(C)C(=O)O)cc1
...
```

### 使用示例

```python
from src.rnn_workflow import RNNMoleculeGenerator

# 初始化
config = {
    'num_molecules': 500,
    'temperature': 0.7,
    'epochs': 100
}
generator = RNNMoleculeGenerator(config)

# 训练
generator.train('data/training_molecules.csv')

# 生成
molecules = generator.generate(500)

# 保存
generator.save_model('models/rnn_model.keras')
```

### 关键参数说明

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `num_molecules` | int | 500 | 生成的分子数量 |
| `temperature` | float | 0.7 | 采样温度，越高越多样 |
| `top_k` | int | 50 | Top-k采样参数 |
| `epochs` | int | 100 | 训练轮数 |
| `batch_size` | int | 32 | 批大小 |

---

## QSAR活性预测模块

### 文件位置
`src/qsar_engine.py`

### 功能描述
基于机器学习的定量构效关系（QSAR）模型，预测分子对靶点的活性（IC50）。

### 核心类/函数

#### `QSARModel`

**主要方法**:

1. `train(train_data_path, model_path)`
   - 功能: 训练QSAR模型
   - 输入: 训练数据路径、模型保存路径
   - 输出: 训练好的模型

2. `predict(molecules_path, output_path)`
   - 功能: 预测分子活性
   - 输入: 分子数据路径
   - 输出: 预测结果CSV

3. `filter_by_ic50(df, threshold)`
   - 功能: 按IC50筛选分子
   - 输入: 数据框、阈值
   - 输出: 筛选后的数据框

### 输入格式

**训练数据** (`data/ChEMBL_Cleaned_akt1-activities_zip.csv`):
```csv
SMILES,IC50
CCOc1ccc(CC(=O)O)cc1,45.2
CC(C)Cc1ccc(C(C)C(=O)O)cc1,82.1
...
```

**待预测分子** (`results/generated_smiles.csv`):
```csv
SMILES
CCOc1ccc(CC(=O)O)cc1
...
```

### 输出格式

**预测结果** (`results/QSAR/predictions.csv`):
```csv
SMILES,pred_IC50,pIC50,MW,LogP,TPSA
CCOc1ccc(CC(=O)O)cc1,45.23,7.34,180.2,2.3,37.3
...
```

### 使用示例

```python
from src.qsar_engine import QSARModel

# 初始化
qsar = QSARModel()

# 训练
qsar.train(
    train_data_path='data/qsar_training.csv',
    model_path='models/qsar_model.pkl'
)

# 预测
predictions = qsar.predict(
    molecules_path='results/generated.csv',
    output_path='results/qsar_predictions.csv'
)

# 筛选
filtered = qsar.filter_by_ic50(predictions, threshold=100)
```

### 分子描述符

模型使用以下分子描述符：

| 描述符 | 说明 |
|--------|------|
| MW | 分子量 |
| LogP | 脂水分配系数 |
| TPSA | 拓扑极性表面积 |
| NumHDonors | 氢键供体数 |
| NumHAcceptors | 氢键受体数 |
| NumRotatableBonds | 可旋转键数 |
| NumAromaticRings | 芳香环数 |
| FractionCSP3 | sp3碳比例 |

---

## ADMET过滤模块

### 文件位置
`src/admet_engine.py`

### 功能描述
评估分子的吸收、分布、代谢、排泄和毒性（ADMET）性质，筛选具有良好药代动力学性质的分子。

### 核心函数

#### `run_admet_filter(df)`

**参数**:
- `df`: 包含SMILES的数据框

**返回**:
- 通过ADMET筛选的数据框

### ADMET指标

| 指标 | 说明 | 筛选标准 |
|------|------|----------|
| P-gp_substrate | P-gp底物 | 0 (非底物) |
| P-gp_I_inhibitor | P-gp I型抑制剂 | 0 (非抑制剂) |
| P-gp_II_inhibitor | P-gp II型抑制剂 | 0 (非抑制剂) |
| CYP2D6_substrate | CYP2D6底物 | 0 (非底物) |
| CYP3A4_substrate | CYP3A4底物 | 0 (非底物) |
| CYP1A2_inhibitor | CYP1A2抑制剂 | 0 (非抑制剂) |
| CYP2C19_inhibitor | CYP2C19抑制剂 | 0 (非抑制剂) |
| CYP2C9_inhibitor | CYP2C9抑制剂 | 0 (非抑制剂) |
| AMES_toxicity | AMES毒性 | 0 (无毒性) |
| hERG_II_inhibitor | hERG II型抑制剂 | 0 (非抑制剂) |
| Hepatotoxicity | 肝毒性 | 0 (无毒性) |
| Skin_Sensitisation | 皮肤过敏 | 0 (无过敏) |

### 使用示例

```python
from src.admet_engine import run_admet_filter
import pandas as pd

# 读取QSAR结果
qsar_df = pd.read_csv('results/QSAR/predictions.csv')

# ADMET过滤
admet_passed = run_admet_filter(qsar_df)

# 保存结果
admet_passed.to_csv('results/ADMET/passed.csv', index=False)
```

---

## 分子对接模块

### 文件位置
`src/submit_job_docking.py`

### 功能描述
将分子与靶蛋白进行分子对接，预测结合模式和结合能。

### 支持的对接引擎

1. **AutoDock Vina**
   - 经典对接引擎
   - 快速、稳定
   - 输出结合能 (kcal/mol)

2. **GNINA**
   - 基于深度学习的对接
   - 考虑蛋白质柔性
   - 输出CNN评分

3. **DiffDock**
   - 扩散模型对接
   - 高精度
   - 输出RMSD预测

### 核心函数

#### `submit_docking_job(engines, receptor, ligands)`

**参数**:
- `engines`: 对接引擎列表
- `receptor`: 受体PDB文件路径
- `ligands`: 配体文件列表

**返回**:
- 对接结果数据框

### 输入格式

**受体文件**: PDB格式
```
ATOM      1  N   MET A   1      11.281  86.699  94.383  1.00 54.62           N
ATOM      2  CA  MET A   1      10.447  87.839  94.727  1.00 54.62           C
...
```

**配体文件**: SDF或PDBQT格式

### 输出格式

**对接结果** (`results/docking/docking_result.csv`):
```csv
SMILES,binding_affinity,rmsd_lb,rmsd_ub,engine
CCOc1ccc(CC(=O)O)cc1,-8.5,0.0,1.2,vina
...
```

### 使用示例

```python
from src.submit_job_docking import submit_docking_job

# 提交对接任务
results = submit_docking_job(
    engines=['vina', 'gnina'],
    receptor='pdb_files_AKT1/4gv1.pdb',
    ligands=['ligand1.sdf', 'ligand2.sdf']
)

# 保存结果
results.to_csv('results/docking/results.csv', index=False)
```

---

## SA Score评估模块

### 文件位置
`src/sa_scorer.py`

### 功能描述
计算分子的合成可及性评分（SA Score），评估分子的合成难度。

### 核心函数

#### `calculate_sa_score(smiles)`

**参数**:
- `smiles`: 分子SMILES字符串

**返回**:
- SA Score (1-10)

#### `run_sa_filter(df, threshold)`

**参数**:
- `df`: 包含SMILES的数据框
- `threshold`: SA Score阈值

**返回**:
- 筛选后的数据框

### SA Score解释

| 范围 | 难度 | 说明 |
|------|------|------|
| 1-3 | 容易 | 常见结构，易于合成 |
| 3-6 | 中等 | 需要一定的合成步骤 |
| 6-10 | 困难 | 复杂结构，合成困难 |

### 使用示例

```python
from src.sa_scorer import calculate_sa_score, run_sa_filter

# 计算单个分子的SA Score
score = calculate_sa_score('CCOc1ccc(CC(=O)O)cc1')
print(f"SA Score: {score}")

# 批量筛选
import pandas as pd
df = pd.read_csv('results/docking/results.csv')
filtered = run_sa_filter(df, threshold=6.0)
```

---

## 可视化模块

### 文件位置
`src/visualizer.py`

### 功能描述
生成各种可视化图表，帮助分析工作流结果。

### 核心函数

#### `generate_plots(final_df, qsar_results)`

**参数**:
- `final_df`: 最终结果数据框
- `qsar_results`: QSAR结果数据框

**输出**:
- 可视化图表PNG文件

### 生成的图表

1. **QSAR预测散点图**
   - 实际值 vs 预测值
   - 显示R²和RMSE

2. **化学空间分布图**
   - MW vs LogP
   - 颜色编码SA Score

3. **IC50分布直方图**
   - 预测IC50分布
   - 核密度估计曲线

4. **多维度评估雷达图**
   - 活性、结合力、SA Score等

### 使用示例

```python
from src.visualizer import generate_plots
import pandas as pd

# 读取数据
final_df = pd.read_csv('results/SA/sa_filtered.csv')
qsar_results = pd.read_csv('results/QSAR/predictions.csv')

# 生成图表
generate_plots(final_df, qsar_results)
```

---

## 日志模块

### 文件位置
`src/logger.py`

### 功能描述
记录工作流执行过程中的所有操作和状态。

### 核心类

#### `WorkflowLogger`

**主要方法**:

1. `log_step(step_name, status, duration=None, output_path=None)`
   - 记录步骤执行状态

2. `log_error(step_name, error_message)`
   - 记录错误信息

3. `log_summary(total_duration)`
   - 记录执行摘要

### 日志格式

```
2026-03-11 12:00:00 - INFO - [RNN分子生成] - 开始执行
2026-03-11 12:05:30 - INFO - [RNN分子生成] - 执行成功 (耗时: 330.45秒) -> results/generated.csv
```

### 使用示例

```python
from src.logger import WorkflowLogger

# 初始化
logger = WorkflowLogger('logs/workflow.log')

# 记录步骤
logger.log_step("RNN生成", "开始")
# ... 执行操作
logger.log_step("RNN生成", "完成", duration=330.45, output_path="results/generated.csv")

# 记录错误
try:
    # ... 可能出错的操作
    pass
except Exception as e:
    logger.log_error("步骤名称", str(e))
```

---

## 工具模块

### 文件位置
`src/utils.py`

### 功能描述
提供各种实用工具函数。

### 核心函数

#### `init_dirs()`
创建必要的目录结构。

#### `validate_smiles(smiles)`
验证SMILES字符串的有效性。

#### `calculate_molecular_properties(smiles)`
计算分子性质（MW, LogP, TPSA等）。

#### `convert_sdf_to_csv(sdf_path, csv_path)`
将SDF文件转换为CSV格式。

#### `batch_process(items, process_func, batch_size)`
批量处理数据。

### 使用示例

```python
from src.utils import validate_smiles, calculate_molecular_properties

# 验证SMILES
is_valid = validate_smiles('CCOc1ccc(CC(=O)O)cc1')

# 计算分子性质
props = calculate_molecular_properties('CCOc1ccc(CC(=O)O)cc1')
print(f"MW: {props['MW']}, LogP: {props['LogP']}")
```

---

## 模块依赖关系

```
run_workflow.py
    ├── src/logger.py
    ├── src/rnn_workflow.py
    │   └── src/utils.py
    ├── src/qsar_engine.py
    │   └── src/utils.py
    ├── src/admet_engine.py
    │   └── src/utils.py
    ├── src/submit_job_docking.py
    │   └── src/utils.py
    ├── src/sa_scorer.py
    │   └── src/utils.py
    └── src/visualizer.py
        └── src/utils.py
```

---

## API参考

详细的API文档请参考各模块源代码中的docstring注释。

---

## 更新日志

| 版本 | 日期 | 更新内容 |
|------|------|----------|
| 1.0.0 | 2026-03-11 | 初始版本 |

---

## 贡献指南

如需添加新模块或改进现有模块，请参考 `docs/customization_guide.md`。
