import pandas as pd
import numpy as np
import pickle
from rdkit import Chem
from rdkit.Chem import Descriptors

# 加载模型
model_path = './models/ultimate_ensemble_qsar_model.pkl'
with open(model_path, 'rb') as f:
    ensemble_model, selector, imputer, scaler = pickle.load(f)

def smiles_to_mol(smiles):
    """Converts a SMILES string to an RDKit molecule object. Handles invalid SMILES."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            Chem.SanitizeMol(mol)
        return mol
    except:
        return None

def calculate_rdkit_descriptors(mol):
    """Calculates RDKit descriptors for a given molecule."""
    if mol is None:
        return None
    return Descriptors.CalcMolDescriptors(mol)

def predict_ic50(smiles_list):
    """Predicts IC50 values for a list of SMILES strings."""
    # 计算描述符
    descriptors_list = []
    valid_smiles = []
    for smiles in smiles_list:
        mol = smiles_to_mol(smiles)
        if mol is not None:
            desc = calculate_rdkit_descriptors(mol)
            if desc is not None:
                descriptors_list.append(desc)
                valid_smiles.append(smiles)
    
    # 转换为DataFrame
    df_desc = pd.DataFrame(descriptors_list)
    
    # 预处理
    df_desc_imputed = imputer.transform(df_desc)
    df_desc_scaled = scaler.transform(df_desc_imputed)
    df_desc_selected = selector.transform(df_desc_scaled)
    
    # 模型预测
    pIC50_pred = ensemble_model.predict(df_desc_selected)
    
    # 转换为IC50
    IC50_pred = 1e9 * np.power(10, -pIC50_pred)
    
    return valid_smiles, IC50_pred, pIC50_pred

def admet_filter(mol):
    """Performs ADMET filtering on a molecule."""
    if mol is None:
        return False
    
    # 简单的ADMET过滤规则
    # 分子量在200-500之间
    mw = Descriptors.MolWt(mol)
    if mw < 200 or mw > 500:
        return False
    
    # LogP在-2到5之间
    logp = Descriptors.MolLogP(mol)
    if logp < -2 or logp > 5:
        return False
    
    # 氢键供体≤5
    hbd = Descriptors.NumHDonors(mol)
    if hbd > 5:
        return False
    
    # 氢键受体≤10
    hba = Descriptors.NumHAcceptors(mol)
    if hba > 10:
        return False
    
    # 可旋转键≤10
    rb = Descriptors.NumRotatableBonds(mol)
    if rb > 10:
        return False
    
    return True

def main():
    # 读取生成的分子
    all_molecules = pd.read_csv('./results/all_generated_molecules.csv')
    top_100_molecules = pd.read_csv('./results/top_100_generated_molecules.csv')
    
    # 预测all_generated_molecules
    print("Predicting all generated molecules...")
    all_smiles = all_molecules['SMILES'].tolist()
    valid_smiles_all, ic50_all, pic50_all = predict_ic50(all_smiles)
    
    # 创建结果DataFrame
    all_results = pd.DataFrame({
        'SMILES': valid_smiles_all,
        'IC50': ic50_all,
        'pIC50': pic50_all
    })
    
    # ADMET过滤
    print("Performing ADMET filtering...")
    all_results['ADMET_Passed'] = all_results['SMILES'].apply(
        lambda x: admet_filter(smiles_to_mol(x))
    )
    
    # 保存结果
    all_results.to_csv('./results/all_generated_molecules_with_predictions.csv', index=False)
    print("All molecules predictions saved to ./results/all_generated_molecules_with_predictions.csv")
    
    # 预测top_100_generated_molecules
    print("Predicting top 100 generated molecules...")
    top_smiles = top_100_molecules['SMILES'].tolist()
    valid_smiles_top, ic50_top, pic50_top = predict_ic50(top_smiles)
    
    # 创建结果DataFrame
    top_results = pd.DataFrame({
        'SMILES': valid_smiles_top,
        'IC50': ic50_top,
        'pIC50': pic50_top
    })
    
    # ADMET过滤
    top_results['ADMET_Passed'] = top_results['SMILES'].apply(
        lambda x: admet_filter(smiles_to_mol(x))
    )
    
    # 保存结果
    top_results.to_csv('./results/top_100_generated_molecules_with_predictions.csv', index=False)
    print("Top 100 molecules predictions saved to ./results/top_100_generated_molecules_with_predictions.csv")
    
    # 统计通过ADMET过滤的比例
    all_passed = all_results['ADMET_Passed'].sum()
    all_total = len(all_results)
    all_ratio = all_passed / all_total * 100
    
    top_passed = top_results['ADMET_Passed'].sum()
    top_total = len(top_results)
    top_ratio = top_passed / top_total * 100
    
    print(f"\nADMET Filter Results:")
    print(f"All molecules: {all_passed}/{all_total} ({all_ratio:.2f}%) passed")
    print(f"Top 100 molecules: {top_passed}/{top_total} ({top_ratio:.2f}%) passed")

if __name__ == "__main__":
    main()
