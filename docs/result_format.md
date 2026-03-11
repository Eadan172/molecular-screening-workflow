# 结果文件格式说明

本文档详细说明工作流各阶段输出的文件格式和内容。

## 目录结构

```
results/
├── QSAR/
│   ├── top_100_generated_molecules_with_predictions.csv
│   ├── all_generated_molecules_with_predictions.csv
│   └── qsar_training_data.csv
├── ADMET/
│   └── top_100_generated_molecules_only_admet_passed.csv
├── docking/
│   ├── docking_result.csv
│   └── docking_poses/
├── SA/
│   └── sa_filtered.csv
└── evaluation_summary.png
```

## 文件格式详解

### 1. RNN生成结果

**文件**: `results/generated_smiles_from_1994mols.csv`

**格式**:
```csv
SMILES
CCOc1ccc(CC(=O)O)cc1
CC(C)Cc1ccc(C(C)C(=O)O)cc1
...
```

**字段说明**:
- `SMILES`: 生成的分子SMILES表示

**数据量**: 默认生成500个分子

---

### 2. QSAR预测结果

**文件**: `results/QSAR/top_100_generated_molecules_with_predictions.csv`

**格式**:
```csv
SMILES,pred_IC50,pIC50,MW,LogP,TPSA
CCOc1ccc(CC(=O)O)cc1,45.23,7.34,180.2,2.3,37.3
CC(C)Cc1ccc(C(C)C(=O)O)cc1,82.15,7.09,194.2,3.1,37.3
...
```

**字段说明**:
- `SMILES`: 分子SMILES表示
- `pred_IC50`: 预测的IC50值 (nM)
- `pIC50`: -log10(IC50)，值越大活性越好
- `MW`: 分子量 (Molecular Weight)
- `LogP`: 脂水分配系数
- `TPSA`: 拓扑极性表面积

**筛选标准**: IC50 < 100 nM

---

### 3. ADMET过滤结果

**文件**: `results/ADMET/top_100_generated_molecules_only_admet_passed.csv`

**格式**:
```csv
SMILES,pred_IC50,P-gp_substrate,P-gp_I_inhibitor,...,Skin_Sensitisation
CCOc1ccc(CC(=O)O)cc1,45.23,0,0,...,0
...
```

**字段说明**:
- `SMILES`: 分子SMILES表示
- `pred_IC50`: QSAR预测的IC50值
- `P-gp_substrate`: P-gp底物 (0=否, 1=是)
- `P-gp_I_inhibitor`: P-gp I型抑制剂 (0=否, 1=是)
- `P-gp_II_inhibitor`: P-gp II型抑制剂 (0=否, 1=是)
- `CYP2D6_substrate`: CYP2D6底物 (0=否, 1=是)
- `CYP3A4_substrate`: CYP3A4底物 (0=否, 1=是)
- `CYP1A2_inhibitor`: CYP1A2抑制剂 (0=否, 1=是)
- `CYP2C19_inhibitor`: CYP2C19抑制剂 (0=否, 1=是)
- `CYP2C9_inhibitor`: CYP2C9抑制剂 (0=否, 1=是)
- `AMES_toxicity`: AMES毒性 (0=无, 1=有)
- `hERG_II_inhibitor`: hERG II型抑制剂 (0=否, 1=是)
- `Hepatotoxicity`: 肝毒性 (0=无, 1=有)
- `Skin_Sensitisation`: 皮肤过敏 (0=无, 1=有)

**筛选标准**: 所有ADMET指标为0（无问题）

---

### 4. 分子对接结果

**文件**: `results/docking/docking_result.csv`

**格式**:
```csv
SMILES,binding_affinity,rmsd_lb,rmsd_ub
CCOc1ccc(CC(=O)O)cc1,-8.5,0.0,1.2
CC(C)Cc1ccc(C(C)C(=O)O)cc1,-7.8,0.0,2.1
...
```

**字段说明**:
- `SMILES`: 分子SMILES表示
- `binding_affinity`: 结合能 (kcal/mol)，越负越好
- `rmsd_lb`: RMSD下界
- `rmsd_ub`: RMSD上界

**对接构象文件**: `results/docking/docking_poses/`
- 格式: PDB或PDBQT
- 命名: `{ligand_name}_docked.pdb`

**筛选标准**: binding_affinity < -7.0 kcal/mol

---

### 5. SA Score评估结果

**文件**: `results/SA/sa_filtered.csv`

**格式**:
```csv
SMILES,sa_score,binding_affinity,pred_IC50,rank
CCOc1ccc(CC(=O)O)cc1,3.2,-8.5,45.23,1
CC(C)Cc1ccc(C(C)C(=O)O)cc1,4.1,-7.8,82.15,2
...
```

**字段说明**:
- `SMILES`: 分子SMILES表示
- `sa_score`: 合成可及性评分 (1-10)
  - 1-3: 易于合成
  - 3-6: 中等难度
  - 6-10: 难以合成
- `binding_affinity`: 对接结合能
- `pred_IC50`: QSAR预测IC50
- `rank`: 综合排名

**排序标准**: 按SA Score升序排列

---

### 6. 可视化报告

**文件**: `results/evaluation_summary.png`

**包含图表**:

1. **QSAR预测散点图**
   - X轴: 实际IC50
   - Y轴: 预测IC50
   - 显示R²值和回归线

2. **化学空间分布图**
   - X轴: 分子量 (MW)
   - Y轴: LogP
   - 颜色: SA Score

3. **IC50分布直方图**
   - X轴: IC50值
   - Y轴: 分子数量
   - 包含核密度估计曲线

4. **SA Score箱线图**
   - 显示最终候选分子的SA Score分布

---

## 日志文件格式

**文件**: `logs/workflow_YYYYMMDD_HHMMSS.log`

**格式示例**:
```
2026-03-11 12:00:00 - INFO - [初始化] - 目录结构创建完成
2026-03-11 12:00:05 - INFO - [RNN分子生成] - 开始执行
2026-03-11 12:05:30 - INFO - [RNN分子生成] - 执行成功 (耗时: 325.45秒) -> results/generated_smiles.csv
2026-03-11 12:05:35 - INFO - [QSAR活性预测] - 开始执行
2026-03-11 12:06:10 - INFO - [QSAR活性预测] - 执行成功 (耗时: 35.23秒) -> results/QSAR/top_100.csv
...
2026-03-11 12:15:00 - INFO - 工作流执行完成！总耗时: 900.12秒
```

**日志级别**:
- `INFO`: 正常执行信息
- `WARNING`: 警告信息（不影响执行）
- `ERROR`: 错误信息（影响执行）

---

## 中间文件

工作流会保留中间结果文件，便于调试和恢复：

```
results/
├── all_generated_molecules.csv              # RNN生成的所有分子
├── all_generated_molecules_with_predictions.csv  # QSAR预测的所有结果
├── train_cleaned.csv                        # 清洗后的训练数据
└── qsar_validation_data.csv                 # QSAR验证数据
```

---

## 数据统计

每个阶段的数据量变化：

| 阶段 | 输入数量 | 输出数量 | 筛选率 |
|------|---------|---------|--------|
| RNN生成 | - | 500 | - |
| QSAR筛选 | 500 | 100 | 80% |
| ADMET过滤 | 100 | ~50 | ~50% |
| 分子对接 | ~50 | ~30 | ~40% |
| SA Score | ~30 | ~20 | ~33% |

---

## 结果解读指南

### 1. 如何判断分子质量？

**综合评估指标**:
- IC50 < 50 nM: 优秀活性
- binding_affinity < -8.0 kcal/mol: 强结合力
- SA Score < 4.0: 易于合成
- ADMET全通过: 良好的药代动力学性质

### 2. 如何选择候选分子？

**推荐排序**:
1. 首先看IC50预测值（活性）
2. 其次看binding_affinity（结合力）
3. 再看SA Score（合成难度）
4. 最后看ADMET性质（安全性）

### 3. 如何理解可视化图表？

**化学空间分布图**:
- 聚集在一起的分子具有相似性质
- 分布广泛的分子库具有更好的多样性

**IC50分布**:
- 左偏分布表示大部分分子活性较好
- 右偏分布表示需要进一步优化

---

## 导出格式

结果文件支持多种导出格式：

1. **CSV**: 默认格式，可用Excel打开
2. **SDF**: 分子结构文件，可用ChemDraw打开
3. **PDB**: 3D结构文件，可用PyMOL打开
4. **JSON**: 结构化数据，便于程序处理

**转换命令**:
```bash
# CSV转SDF
python scripts/convert_csv_to_sdf.py results/SA/sa_filtered.csv

# CSV转JSON
python -c "import pandas as pd; pd.read_csv('results.csv').to_json('results.json')"
```

---

## 注意事项

1. **文件编码**: 所有CSV文件使用UTF-8编码
2. **数值精度**: IC50保留2位小数，其他保留1位小数
3. **缺失值**: 使用空字符串或NaN表示
4. **文件大小**: 单个CSV文件建议不超过10MB

---

## 联系支持

如有疑问，请：
1. 查看日志文件 `logs/workflow_*.log`
2. 参考文档 `docs/`
3. 提交Issue: https://github.com/yourusername/molecular-screening-workflow/issues
