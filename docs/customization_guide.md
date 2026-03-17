# 定制化与扩展指南

本文档说明如何定制和扩展分子筛选工作流。

## 目录

1. [配置文件定制](#配置文件定制)
2. [添加新的筛选标准](#添加新的筛选标准)
3. [集成新的对接引擎](#集成新的对接引擎)
4. [自定义机器学习模型](#自定义机器学习模型)
5. [扩展工作流步骤](#扩展工作流步骤)
6. [性能优化](#性能优化)

---

## 配置文件定制

### 修改配置文件

编辑 `config/workflow_config.yaml`:

```yaml
# 示例：修改QSAR筛选标准
qsar:
  ic50_threshold: 50  # 从100改为50，更严格
  top_n: 200          # 从100改为200，保留更多分子

# 示例：修改RNN生成参数
rnn:
  num_molecules: 1000  # 生成更多分子
  temperature: 0.8     # 提高多样性
  top_k: 100           # 增加采样范围
```

### 创建多个配置文件

针对不同靶点或场景创建专用配置：

```bash
config/
├── workflow_config.yaml          # 默认配置
├── akt1_config.yaml              # AKT1靶点专用
├── egfr_config.yaml              # EGFR靶点专用
└── high_throughput_config.yaml   # 高通量筛选配置
```

使用自定义配置：
```bash
python run_workflow.py --config config/akt1_config.yaml
```

---

## 添加新的筛选标准

### 示例：添加Lipinski五规则筛选

1. 创建新模块 `src/lipinski_filter.py`:

```python
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

def lipinski_filter(df):
    """应用Lipinski五规则筛选"""
    def check_lipinski(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        
        # Lipinski五规则
        return (mw <= 500 and 
                logp <= 5 and 
                hbd <= 5 and 
                hba <= 10)
    
    df['lipinski_pass'] = df['SMILES'].apply(check_lipinski)
    filtered_df = df[df['lipinski_pass']].copy()
    
    return filtered_df

if __name__ == "__main__":
    df = pd.read_csv('input.csv')
    filtered = lipinski_filter(df)
    filtered.to_csv('lipinski_filtered.csv', index=False)
```

2. 在 `run_workflow.py` 中集成:

```python
from src.lipinski_filter import lipinski_filter

def run_lipinski_filter(self):
    """步骤X: Lipinski五规则筛选"""
    self.logger.log_step("Lipinski筛选", "开始")
    
    start_time = time.time()
    
    try:
        input_df = self.results.get('admet_df')
        filtered_df = lipinski_filter(input_df)
        
        duration = time.time() - start_time
        output_path = 'results/lipinski_filtered.csv'
        filtered_df.to_csv(output_path, index=False)
        
        self.logger.log_step("Lipinski筛选", "完成", duration, output_path)
        return filtered_df
    except Exception as e:
        self.logger.log_step("Lipinski筛选", f"失败: {str(e)}")
        raise
```

3. 更新配置文件:

```yaml
# 添加Lipinski配置
lipinski:
  enabled: true
  mw_threshold: 500
  logp_threshold: 5
  hbd_threshold: 5
  hba_threshold: 10
```

---

## 集成新的对接引擎

### 示例：集成Gold对接软件

1. 创建对接引擎接口 `src/gold_docking.py`:

```python
import subprocess
import os

class GoldDocking:
    def __init__(self, config):
        self.gold_path = config.get('gold_path', 'gold')
        self.config_file = config.get('gold_config')
    
    def prepare_ligand(self, smiles):
        """准备配体文件"""
        # 使用RDKit转换SMILES为MOL2
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        
        mol2_path = f"ligand_{hash(smiles)}.mol2"
        Chem.MolToMol2File(mol, mol2_path)
        return mol2_path
    
    def run_docking(self, receptor, ligand, output_dir):
        """运行Gold对接"""
        cmd = [
            self.gold_path,
            '-c', self.config_file,
            '-r', receptor,
            '-l', ligand,
            '-o', output_dir
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        return self.parse_results(output_dir)
    
    def parse_results(self, output_dir):
        """解析Gold对接结果"""
        # 解析Gold输出文件
        results = []
        # ... 解析逻辑
        return results
```

2. 在 `run_workflow.py` 中添加选项:

```python
parser.add_argument(
    '--docking-engine',
    choices=['vina', 'gnina', 'diffdock', 'gold'],
    default='vina',
    help='选择对接引擎'
)
```

3. 更新配置文件:

```yaml
docking:
  engines:
    - gold
  
  gold:
    path: /path/to/gold
    config_file: config/gold_config.conf
    genetic_algorithm:
      population_size: 100
      max_operations: 100000
```

---

## 自定义机器学习模型

### 示例：使用XGBoost替代Random Forest

1. 创建新模型 `src/qsar_xgboost.py`:

```python
import xgboost as xgb
from sklearn.model_selection import cross_val_score
import pandas as pd
import numpy as np

class XGBoostQSAR:
    def __init__(self, config):
        self.model = xgb.XGBRegressor(
            n_estimators=config.get('n_estimators', 100),
            max_depth=config.get('max_depth', 6),
            learning_rate=config.get('learning_rate', 0.1),
            random_state=42
        )
    
    def train(self, X, y):
        """训练XGBoost模型"""
        self.model.fit(X, y)
        
        # 交叉验证
        scores = cross_val_score(self.model, X, y, cv=5, scoring='r2')
        print(f"交叉验证R²: {scores.mean():.4f} (+/- {scores.std():.4f})")
        
        return self.model
    
    def predict(self, X):
        """预测活性"""
        return self.model.predict(X)
    
    def save(self, path):
        """保存模型"""
        self.model.save_model(path)
    
    def load(self, path):
        """加载模型"""
        self.model.load_model(path)

def train_xgboost_qsar(train_data_path, model_path):
    """训练XGBoost QSAR模型"""
    df = pd.read_csv(train_data_path)
    
    # 准备特征和标签
    X = df.drop(['SMILES', 'IC50'], axis=1)
    y = df['IC50']
    
    # 训练模型
    model = XGBoostQSAR({'n_estimators': 200})
    model.train(X, y)
    model.save(model_path)
    
    return model
```

2. 更新配置文件:

```yaml
qsar:
  model_type: xgboost  # 或 random_forest
  
  xgboost:
    n_estimators: 200
    max_depth: 8
    learning_rate: 0.05
```

---

## 扩展工作流步骤

### 添加新的工作流步骤

1. 在 `run_workflow.py` 中添加新方法:

```python
def run_custom_analysis(self):
    """步骤X: 自定义分析"""
    self.logger.log_step("自定义分析", "开始")
    
    start_time = time.time()
    
    try:
        # 你的自定义逻辑
        input_data = self.results.get('previous_step_output')
        result = self.custom_analysis_function(input_data)
        
        duration = time.time() - start_time
        output_path = 'results/custom_analysis.csv'
        result.to_csv(output_path, index=False)
        
        self.logger.log_step("自定义分析", "完成", duration, output_path)
        return result
    except Exception as e:
        self.logger.log_step("自定义分析", f"失败: {str(e)}")
        raise
```

2. 在 `run_full_workflow` 中调用:

```python
def run_full_workflow(self):
    # ... 现有步骤
    
    # 添加新步骤
    custom_result = self.run_custom_analysis()
    self.results['custom_result'] = custom_result
    
    # ... 后续步骤
```

---

## 性能优化

### 1. 并行处理

在配置文件中启用并行:

```yaml
parallel:
  enabled: true
  num_workers: 4  # 根据CPU核心数调整
```

在代码中使用:

```python
from joblib import Parallel, delayed

def parallel_docking(ligands, receptor, n_jobs=4):
    """并行对接"""
    results = Parallel(n_jobs=n_jobs)(
        delayed(dock_single_ligand)(ligand, receptor)
        for ligand in ligands
    )
    return results
```

### 2. 内存优化

使用生成器处理大数据:

```python
def batch_generator(data, batch_size=100):
    """批量数据生成器"""
    for i in range(0, len(data), batch_size):
        yield data[i:i+batch_size]

# 使用
for batch in batch_generator(large_dataset):
    process_batch(batch)
```

### 3. 缓存机制

缓存计算结果:

```python
from functools import lru_cache

@lru_cache(maxsize=1000)
def calculate_descriptors(smiles):
    """缓存分子描述符计算"""
    mol = Chem.MolFromSmiles(smiles)
    return Descriptors.CalcMolDescriptors(mol)
```

---

## 添加新的可视化

### 自定义图表

创建 `src/custom_visualizer.py`:

```python
import matplotlib.pyplot as plt
import seaborn as sns

def plot_custom_analysis(df, output_path):
    """自定义可视化"""
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # 图1: 活性分布
    sns.histplot(df['pred_IC50'], ax=axes[0, 0])
    axes[0, 0].set_title('IC50 Distribution')
    
    # 图2: 性质相关性
    sns.heatmap(df.corr(), ax=axes[0, 1], cmap='coolwarm')
    axes[0, 1].set_title('Property Correlation')
    
    # 图3: 化学空间
    sns.scatterplot(data=df, x='MW', y='LogP', 
                    hue='sa_score', ax=axes[1, 0])
    axes[1, 0].set_title('Chemical Space')
    
    # 图4: 排名分布
    df_sorted = df.sort_values('rank').head(20)
    sns.barplot(data=df_sorted, x='rank', y='pred_IC50', ax=axes[1, 1])
    axes[1, 1].set_title('Top 20 Candidates')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
```

---

## 测试与验证

### 单元测试

创建 `tests/test_custom_module.py`:

```python
import pytest
from src.custom_module import custom_function

def test_custom_function():
    """测试自定义函数"""
    input_data = "test_smiles"
    result = custom_function(input_data)
    
    assert result is not None
    assert len(result) > 0

if __name__ == "__main__":
    pytest.main([__file__, '-v'])
```

运行测试:
```bash
python -m pytest tests/ -v
```

---

## 最佳实践

1. **版本控制**: 为自定义模块创建独立分支
2. **文档更新**: 更新README和docs文档
3. **配置分离**: 使用独立的配置文件
4. **错误处理**: 添加详细的异常处理和日志
5. **性能测试**: 对新功能进行性能基准测试

---

## 获取帮助

- 查看现有模块源代码
- 参考API文档
- 提交Issue讨论
- 查看示例代码 `examples/`

---

**祝你定制愉快！** 🚀
