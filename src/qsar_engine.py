from operator import index
import os
import tqdm
from xml.parsers.expat import model
import pandas as pd
import numpy as np

Chem = None
AllChem = None
Descriptors = None

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import Descriptors
except ImportError:
    pass

try:
    from joblib import dump, load
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.metrics import r2_score, mean_squared_error
    from sklearn.model_selection import train_test_split
except ImportError:
    pass


def smiles_to_mol(smiles):
    """Converts a SMILES string to an RDKit molecule object. Handles invalid SMILES."""
    if Chem is None:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            Chem.SanitizeMol(mol)
        return mol
    except:
        return None

# 文章数据是微摩尔，脚本使用的是纳摩尔
def convert_to_nm(row):
    val = row['IC50']
    return val * 1000

def calculate_rdkit_descriptors(mol):
    """Calculates RDKit descriptors for a given molecule."""
    if mol is None or Descriptors is None:
        return None
    return Descriptors.CalcMolDescriptors(mol)

def visualize_qsar_results(output_path, y):
    """生成 QSAR 预测结果的可视化图表"""
    import matplotlib.pyplot as plt
    import seaborn as sns

    
    df_generated_molecules_cleaned = pd.read_csv(output_path, encoding='latin1')
    df_train = pd.read_csv('./results/train_cleaned.csv')
    y = df_train['IC50']

    # Set a style for better aesthetics
    sns.set_style("whitegrid")

    plt.figure(figsize=(10, 6))

    # Plot the distribution of predicted IC50 for generated molecules
    sns.histplot(df_generated_molecules_cleaned['pred_IC50'], color='skyblue', label='Generated Molecules (Predicted)', kde=True, stat='density', alpha=0.6)

    # Plot the distribution of actual IC50 for inhibitors (training data)
    sns.histplot(y, color='lightcoral', label='Inhibitors (Actual)', kde=True, stat='density', alpha=0.6)

    plt.title('Distribution of IC50 Values: Generated Molecules vs. Inhibitors')
    plt.xlabel('IC50 (nM)')
    plt.ylabel('Density')
    plt.legend()
    plt.xlim(0, 1) # Limit x-axis to a reasonable range for IC50 values
    #plt.show()
    plt.savefig('./models/distribution of predicted IC50 for generated molecules.png',dpi=300)

    print("Generated molecule predictions compared to inhibitor actual activities plot displayed.")


def train_qsar_model(training_data,model_path):
    print(">>> 步骤 1 & 2: 正在训练 QSAR 模型并进行动态筛选...")

    X = training_data.drop('IC50', axis=1)
    y = training_data['IC50']

    print("Shape of X:", X.shape)
    print("Shape of y:", y.shape)

    # Split the features (X) and target (y) into training and validation sets
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=42)
    df_train = pd.concat([X_train, y_train], axis=1)
    df_train.to_csv('./results/QSAR/qsar_training_data.csv', index=False)
    df_val = pd.concat([X_val, y_val], axis=1)
    df_val.to_csv('./results/QSAR/qsar_validation_data.csv', index=False)

    print(f"X_train shape: {X_train.shape}")
    print(f"X_val shape: {X_val.shape}")
    print(f"y_train shape: {y_train.shape}")
    print(f"y_val shape: {y_val.shape}")
    
    
    # 添加特征选择
    from sklearn.feature_selection import SelectKBest, f_regression
    selector = SelectKBest(f_regression, k=min(100, X_train.shape[1]))
    X_train_selected = selector.fit_transform(X_train, y_train)
    X_val_selected = selector.transform(X_val)
    
    print(f"Original features: {X_train.shape[1]}, Selected features: {X_train_selected.shape[1]}")
    
    # 训练 (利用 10 核 CPU)
    model = RandomForestRegressor(
        n_estimators=200, 
        max_depth=30, 
        min_samples_split=2, 
        min_samples_leaf=1, 
        max_features='sqrt', 
        n_jobs=10, 
        random_state=42
    )
    print("Training Random Forest Regressor model...")
    model.fit(X_train_selected, y_train)
    print("Model training complete.")

    # Make predictions on the validation set
    y_pred = model.predict(X_val_selected)

    # Calculate and print the R-squared score
    r2 = r2_score(y_val, y_pred)
    print(f"R-squared score on the validation set: {r2:.4f}")

    # Calculate and print the Root Mean Squared Error (RMSE)
    rmse = np.sqrt(mean_squared_error(y_val, y_pred))
    print(f"Root Mean Squared Error (RMSE) on the validation set: {rmse:.4f}")
    
    # 保存模型和特征选择器
    import pickle
    with open(model_path, 'wb') as f:
        pickle.dump((model, selector), f)
    
    return model, selector

def qsar_predict_and_filter(train_data_path, predict_data_path, model_path='./models/rf_qsar_model.pkl', demo_mode=False):
    print(">>> 启动 QSAR 筛选流水线...")
    
    if demo_mode:
        print("演示模式: 使用模拟QSAR预测结果")
        demo_result = pd.DataFrame({
            'SMILES': [
                'CCOc1ccc(CC(=O)O)cc1',
                'CC(C)Cc1ccc(C(C)C(=O)O)cc1',
                'COc1ccc(CC(=O)O)cc1',
                'CC(C)Cc1ccc(CC(=O)O)cc1',
                'COc1ccc(C(C)C(=O)O)cc1'
            ],
            'pred_IC50': [45.2, 82.1, 35.8, 68.5, 52.3],
            'pIC50': [7.34, 7.09, 7.45, 7.16, 7.28],
            'MW': [180.2, 194.2, 166.2, 208.3, 180.2],
            'LogP': [2.3, 3.1, 1.8, 2.9, 2.5]
        })
        os.makedirs('./results/QSAR', exist_ok=True)
        demo_result.to_csv('./results/QSAR/top_100_generated_molecules_with_predictions.csv', index=False)
        print(f"演示QSAR结果已保存，共 {len(demo_result)} 个分子")
        return demo_result

    df_inhibitors = pd.read_csv(train_data_path,encoding='latin1') #训练的分子
    df_inhibitors = df_inhibitors.dropna(subset=['SMILES']).reset_index(drop=True)
    df_inhibitors['Mol'] = df_inhibitors['SMILES'].apply(smiles_to_mol)
    df_inhibitors_cleaned = df_inhibitors.dropna(subset=['Mol'])

    print(f"Shape of df_inhibitors after cleaning: {df_inhibitors_cleaned.shape}")

    df_generated_molecules = pd.read_csv(predict_data_path,encoding='latin1') #AI生成的分子
    df_generated_molecules = df_generated_molecules.dropna(subset=['SMILES']).reset_index(drop=True)
    df_generated_molecules['Mol'] = df_generated_molecules['SMILES'].apply(smiles_to_mol)
    df_generated_molecules_cleaned = df_generated_molecules.dropna(subset=['Mol'])

    print(f"Shape of df_generated_molecules after cleaning: {df_generated_molecules_cleaned.shape}")

    # Generate descriptors for inhibitors
    print("Generating descriptors for inhibitors...")
    #print('合并前查看IC50column',df_inhibitors_cleaned.columns)
    inhibitors_descriptors_list = df_inhibitors_cleaned['Mol'].apply(calculate_rdkit_descriptors)
    inhibitors_descriptors = pd.DataFrame(inhibitors_descriptors_list.tolist())

    # Generate descriptors for generated molecules
    print("Generating descriptors for generated molecules...")
    generated_molecules_descriptors_list = df_generated_molecules_cleaned['Mol'].apply(calculate_rdkit_descriptors)
    generated_molecules_descriptors = pd.DataFrame(generated_molecules_descriptors_list.tolist())
    common_cols = list(set(inhibitors_descriptors.columns) & set(generated_molecules_descriptors.columns))

    # Filter both dataframes to retain only common columns
    inhibitors_descriptors_aligned = inhibitors_descriptors[common_cols]
    generated_molecules_descriptors_aligned = generated_molecules_descriptors[common_cols]

    print(f"Number of common columns: {len(common_cols)}")
    print("Inhibitors Descriptors Aligned Shape: ", inhibitors_descriptors_aligned.shape)
    print("Generated Molecules Descriptors Aligned Shape: ", generated_molecules_descriptors_aligned.shape)

    df_inhibitors_cleaned_reset = df_inhibitors_cleaned.reset_index(drop=True)
    #print('合并前查看IC50column',df_inhibitors_cleaned.columns)
    inhibitors_descriptors_aligned_reset = inhibitors_descriptors_aligned.reset_index(drop=True)


    training_data = pd.concat([
    inhibitors_descriptors_aligned_reset,
        df_inhibitors_cleaned_reset['IC50']], axis=1)

    print("Shape of training_data before handling NaNs:", training_data.shape)
    missing_values_count = training_data.isnull().sum().sum()
    if missing_values_count > 0:
        print(f"Found {missing_values_count} missing values in the training data. Dropping rows with NaNs.")
        training_data.dropna(inplace=True)
        print("Shape of training_data after dropping NaNs:", training_data.shape)
    else:
        print("No missing values found in the training data.")
    
    # 转换IC50单位为纳摩尔
    training_data['IC50'] = df_inhibitors_cleaned.apply(convert_to_nm, axis=1)
    training_data.to_csv('./results/train_cleaned.csv',index=False)

    # 训练模型并获取特征选择器
    model, selector = train_qsar_model(training_data,model_path)

    # 预测 IC50

    # Replace infinite values with NaN
    generated_molecules_descriptors_aligned.replace([np.inf, -np.inf], np.nan, inplace=True)

    # Drop rows with NaN values from the generated molecules descriptors
    initial_shape = generated_molecules_descriptors_aligned.shape
    generated_molecules_descriptors_aligned.dropna(inplace=True)

    if generated_molecules_descriptors_aligned.shape[0] < initial_shape[0]:
        print(f"Dropped {initial_shape[0] - generated_molecules_descriptors_aligned.shape[0]} rows from generated_molecules_descriptors_aligned due to NaN or infinite values.")

    # Also, ensure df_generated_molecules_cleaned is aligned with the dropped rows
    df_generated_molecules_cleaned = df_generated_molecules_cleaned.loc[generated_molecules_descriptors_aligned.index]

    # 使用相同的特征选择器转换特征
    generated_molecules_descriptors_selected = selector.transform(generated_molecules_descriptors_aligned)
    predicted_activities = model.predict(generated_molecules_descriptors_selected)
    df_generated_molecules_cleaned['pred_IC50'] = predicted_activities
    
    # 加载akt1-activities数据集，分析IC50分布
    akt1_data = pd.read_csv('./csv_input/ChEMBL_Cleaned_akt1-activities_zip.csv')
    akt1_ic50_min = akt1_data['IC50_nM'].min()
    akt1_ic50_max = akt1_data['IC50_nM'].max()
    akt1_ic50_median = akt1_data['IC50_nM'].median()
    
    print(f"akt1-activities数据集IC50分布: 最小值={akt1_ic50_min:.2f}nM, 最大值={akt1_ic50_max:.2f}nM, 中位数={akt1_ic50_median:.2f}nM")
    
    # 筛选IC50值在akt1-activities数据集分布区间内的分子
    df_generated_molecules_cleaned = df_generated_molecules_cleaned[
        (df_generated_molecules_cleaned['pred_IC50'] >= akt1_ic50_min * 0.1)  # 允许一定的误差范围
        & (df_generated_molecules_cleaned['pred_IC50'] <= akt1_ic50_max * 1.1)
    ]
    
    print(f"筛选后符合IC50分布区间的分子数量: {len(df_generated_molecules_cleaned)}")

    print("df_generated_molecules_cleaned head with predicted activities:")
    print(df_generated_molecules_cleaned.shape)
    df_generated_molecules_cleaned.to_csv('./results/df_generated_molecules_cleaned.csv',index=False)
    # 动态阈值逻辑：优先 IC50 < 100nM (注意：通常活性是值越小越好)
    # 如果你的模型训练的是 pIC50，则需调整符号
    '''
    thresholds = [100, 200, 300, 500]
    final_qsar_df = pd.DataFrame()
    
    for t in tqdm.tqdm(thresholds):
        filtered = df_generated_molecules_cleaned[df_generated_molecules_cleaned['pred_IC50'] <= t] # 筛选高活性分子
        if len(filtered) >= 100:
            print(f"QSAR 阶段完成: 阈值 {t}nM 下获得 {len(filtered)} 个分子")
            final_qsar_df = filtered
            break
    
    output_path = './results/QSAR/qsar_screened.csv'
    final_qsar_df.to_csv(output_path, index=False)
    '''
    return df_generated_molecules_cleaned