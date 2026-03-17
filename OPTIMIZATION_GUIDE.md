# 分子筛选工作流模型优化指南

## 概述

本指南详细说明了如何优化分子筛选工作流中的RNN和QSAR模型，以达到以下性能目标：

- **RNN模型**: accuracy >= 0.95
- **QSAR模型**: R² > 0.8, RMSE < 0.2

## 数据准备

### 数据合并与去重

已完成数据合并，生成了两个合并后的数据集：

1. **RNN训练数据** (`data/merged_rnn_data.csv`):
   - 原始数据: ChEMBL_Cleaned_Finall_test_zipfile.csv (1,032,632条) + Anticancer_compounds@1994_from_paper.csv (1,994条)
   - 合并后: 1,031,551条
   - 去重后: 538,871条
   - 去除重复: 492,680条

2. **QSAR训练数据** (`data/merged_qsar_data.csv`):
   - 原始数据: ChEMBL_Cleaned_akt1-activities_zip.csv (3,713条) + Anticancer_compounds@1994_from_paper.csv (1,994条)
   - 合并后: 5,707条
   - 去重后: 3,466条
   - 去除重复: 2,241条

### 数据预处理说明

- **RNN数据**: 仅保留SMILES列，用于分子生成模型训练
- **QSAR数据**: 
  - 保留SMILES和IC50列
  - 将Anticancer_compounds数据的IC50从微摩尔(μM)转换为纳摩尔(nM)
  - 确保所有IC50值单位统一为nM

## 模型优化详情

### RNN模型优化

#### 原始架构
- 3层LSTM (256, 256, 128)
- 单向LSTM
- Embedding维度: 256
- Dropout: 0.2
- 学习率: 0.001
- Batch size: 64
- Epochs: 50

#### 优化后架构
- **3层双向LSTM (Bi-LSTM)** (512, 512, 256)
- **Layer Normalization** 层
- Embedding维度: 512
- Dropout: 0.1 (更低的dropout以提高准确率)
- 学习率: 0.0005 (更小的学习率)
- Batch size: 128 (更大的batch size)
- Epochs: 200 (更多训练轮数)
- 验证集: 10%用于验证

#### 训练策略
- EarlyStopping: 监控accuracy，patience=10，min_delta=0.001
- ReduceLROnPlateau: 监控accuracy，factor=0.5，patience=5
- ModelCheckpoint: 保存最佳模型

#### 关键改进点
1. **双向LSTM**: 捕获前后文信息，提高序列建模能力
2. **层归一化**: 稳定训练过程，加速收敛
3. **更大的隐藏层**: 512维替代256维，增加模型容量
4. **优化的训练参数**: 更小的学习率、更大的batch size
5. **验证集监控**: 防止过拟合

### QSAR模型优化

#### 原始架构
- 单模型: RandomForestRegressor
- n_estimators: 200
- max_depth: 30
- 特征选择: SelectKBest (k=100)
- 无标准化

#### 优化后架构
- **集成模型**: VotingRegressor (RandomForest + GradientBoosting)
- **Pipeline**: 标准化 -> 特征选择 -> 集成预测
- RandomForest参数: n_estimators=300, max_depth=30
- GradientBoosting参数: n_estimators=300, max_depth=7, learning_rate=0.1
- 特征选择: SelectKBest (k=150)
- **标准化**: StandardScaler
- **对数变换**: 对IC50进行log1p变换

#### 关键改进点
1. **集成学习**: 结合RandomForest和GradientBoosting的优势
2. **数据标准化**: 使特征尺度一致，提高模型性能
3. **对数变换**: IC50通常呈偏态分布，对数变换使其更接近正态分布
4. **更多特征**: 从100个增加到150个特征
5. **更多树**: n_estimators从200增加到300

## 使用方法

### 1. 仅数据合并（已完成）

```bash
python merge_data_only.py
```

### 2. 完整模型优化（推荐）

运行完整的优化流程，包括数据合并、RNN模型训练和QSAR模型训练：

```bash
python optimize_models.py
```

### 3. 使用原始工作流

配置文件已更新为使用合并的数据。可以直接运行原始工作流：

```bash
# 使用更新后的配置
python run_workflow.py
```

## 文件结构

```
molecular-screening-workflow/
├── data/
│   ├── merged_rnn_data.csv          # 合并后的RNN训练数据
│   ├── merged_qsar_data.csv         # 合并后的QSAR训练数据
│   └── ... (原始数据文件)
├── models/
│   ├── selfies_generator_rnn.keras  # 优化后的RNN模型
│   ├── selfies_generator_rnn_performance.json  # RNN性能记录
│   ├── ultimate_ensemble_qsar_model.pkl  # 优化后的QSAR模型
│   └── ultimate_ensemble_qsar_model_performance.json  # QSAR性能记录
├── src/
│   ├── rnn_workflow.py  # 更新后的RNN工作流
│   ├── qsar_engine.py   # 更新后的QSAR引擎
│   └── ...
├── config/
│   └── workflow_config.yaml  # 更新后的配置文件
├── optimize_models.py    # 完整优化脚本
├── merge_data_only.py    # 仅数据合并脚本
└── OPTIMIZATION_GUIDE.md  # 本文档
```

## 性能监控

### RNN模型性能指标
- 训练准确率 (accuracy)
- 验证准确率 (val_accuracy)
- 训练损失 (loss)
- 验证损失 (val_loss)

### QSAR模型性能指标
- R²分数 (决定系数)
- RMSE (均方根误差)
- 交叉验证R²分数

性能记录会自动保存在模型同名的`_performance.json`文件中。

## 预期结果

### RNN模型
- 目标: accuracy >= 0.95
- 预期训练时间: 取决于硬件，建议使用GPU

### QSAR模型
- 目标: R² > 0.8, RMSE < 0.2
- 预期训练时间: 10-30分钟（取决于CPU核心数）

## 注意事项

1. **硬件要求**:
   - RNN训练建议使用GPU（CUDA支持）
   - QSAR训练可以使用CPU，建议多核

2. **依赖库**:
   - 确保已安装所有依赖: `pip install -r requirements.txt`
   - 主要依赖: TensorFlow/Keras, scikit-learn, RDKit, SELFIES

3. **内存要求**:
   - RNN训练可能需要大量内存（建议16GB+）
   - 如遇内存问题，可以减小batch_size

4. **训练时间**:
   - 完整训练可能需要数小时
   - 可以调整epochs参数减少训练时间

## 故障排除

### RNN模型训练问题
- **内存不足**: 减小batch_size或减少LSTM层维度
- **训练不收敛**: 尝试调整学习率或增加dropout
- **准确率过低**: 增加训练轮数或调整模型架构

### QSAR模型训练问题
- **R²过低**: 尝试更多特征或调整集成模型权重
- **RMSE过高**: 检查数据预处理，确保IC50单位正确
- **训练过慢**: 减少n_estimators或使用更少的CPU核心

## 下一步

模型训练完成后，可以：

1. 查看性能记录JSON文件
2. 运行完整的分子筛选工作流
3. 根据实际结果进一步调优超参数
4. 生成新的分子并进行筛选

## 联系方式

如有问题，请查看日志文件或联系开发团队。
