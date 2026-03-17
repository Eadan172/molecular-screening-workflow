import argparse
from importlib.resources import path
import sys
from qsar_engine import qsar_predict_and_filter
from admet_engine import run_admet_filter
from docking_engine import parse_and_save_results, run_vina_parallel, run_gnina_parallel, run_diffdock, prepare_docking_input
from submit_job_docking import submit_docking_job
import pandas as pd
import os

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from joblib import dump, load

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.model_selection import train_test_split


# 1. 目录初始化 
def init_dirs():
    dirs = ['./models', './results/QSAR', './results/ADMET', './results/docking', './results/SA']
    for d in dirs:
        os.makedirs(d, exist_ok=True)

def smiles_to_mol(smiles):
    """Converts a SMILES string to an RDKit molecule object. Handles invalid SMILES."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            Chem.SanitizeMol(mol)
        return mol
    except:
        return None
    
# 5. SA Score 筛选
def run_sa_filter(df):
    from rdkit.Chem import RDConfig
    import sys
    sys.path.append(f'{RDConfig.RDContribDir}/SA_Score')
    import sascorer

    print(">>> 步骤 5: 正在计算 SA Score...")
    df['sa_score'] = df['mol'].apply(lambda x: sascorer.calculateScore(x))
    # 按 SA Score 排序并保存
    sa_df = df.sort_values('sa_score')
    sa_df.to_csv('./results/SA/sa_filtered.csv', index=False)
    return sa_df

# 6. 可视化
def generate_plots(final_df, qsar_results):
    generate_plots(final_df, qsar_results)
    
    plt.figure(figsize=(16, 10))
    
    # 1. QSAR 预测 vs 实际 (散点/回归)
    plt.subplot(2, 2, 1)
    sns.scatterplot(x=y_test, y=y_pred, alpha=0.5)
    plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--')
    plt.title(f"QSAR Training: R²={metrics['r2']:.2f}")
    
    # 2. 化学空间分布 (MW vs LogP)
    plt.subplot(2, 2, 2)
    sns.scatterplot(data=final_df, x='MW', y='LogP', hue='sa_score', palette='viridis')
    plt.title("Chemical Space (MW vs LogP)")

    # 3. IC50 柱状分布
    plt.subplot(2, 2, 3)
    sns.histplot(final_df['pred_IC50'], bins=20, color='skyblue', kde=True)
    plt.title("Predicted IC50 Distribution")

    # 4. ADMET 属性雷达图/箱线图示例 (展示 MW 分布)
    plt.subplot(2, 2, 4)
    sns.boxplot(data=final_df, y='sa_score', color='lightgreen')
    plt.title("SA Score Distribution of Final Hits")

    plt.tight_layout()
    plt.savefig('./results/evaluation_summary.png')
    print("可视化结果已保存至 ./results/evaluation_summary.png")

def main():
    final_qsar_df = pd.read_csv('./results/QSAR/top_100_generated_molecules_with_predictions.csv')
    admet_df = run_admet_filter(final_qsar_df)
    
    parser = argparse.ArgumentParser(description="分阶段分子评价虚筛工作流")
    parser.add_argument('--rnn_train_data', type=str, required=True, help="原始分子库 CSV")
    parser.add_argument('--rnn_generate_data', type=str, required=True, help="需要预测的分子数据 CSV")
    parser.add_argument('--rnn_model_path', type=str, default='./results/models',required=True, help="训练好的 QSAR 模型路径")
    parser.add_argument('--qsar_train_data', type=str, required=True, help="靶点活性数据")
    parser.add_argument('--qsar_generate_data', type=str, required=True, help="需要预测的分子")
    parser.add_argument('--qsar_model_path', type=str, default='./results/models',required=True, help="训练好的 QSAR 模型路径")
    parser.add_argument('--docking_engine', nargs='+', default=['vina'], help="选择对接程序，可选gnina、autodock-vina、diffdock-l")
    args = parser.parse_args()

    qsar_train_data_path = args.train_data
    qsar_predict_data_path = args.predict_data
    model_path = args.model_path
    docking_engine = args.docking_engine

    init_dirs()

    # 1. RNN模型训练及分子生成

    # 2. QSAR
    final_qsar_df = qsar_predict_and_filter(qsar_train_data_path, qsar_predict_data_path, model_path)
    
    # 3. ADMET
    admet_df = run_admet_filter(final_qsar_df)
    
    # 4. docking
    new_job_id = submit_docking_job(docking_engine)
    

    # 5.SA score
    run_sa_filter(pd.read_csv('./results/docking/docking_result.csv'))

    # 6.visualization
    generate_plots(final_df, qsar_results)

if __name__ == "__main__":
    main()