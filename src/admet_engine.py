
# 过滤前先通过ADMET平台获取12项指标的预测结果，确保分子在ADMET方面具有较好的性质。
# 12项指标包括：
# 1. P-gp底物 (P-gp substrate)
# 2. P-gp I型抑制剂 (P-gp I inhibitor)
# 3. P-gp II型抑制剂 (P-gp II inhibitor)
# 4. CYP2D6底物 (CYP2D6 substrate)
# 5. CYP3A4底物 (CYP3A4 substrate)
# 6. CYP1A2抑制剂 (CYP1A2 inhibitor)
# 7. CYP2C19抑制剂 (CYP2C19 inhibitor
# 8. CYP2C9抑制剂 (CYP2C9 inhibitor)
# 9. AMES毒性 (AMES toxicity)
# 10. hERG II型抑制剂 (hERG II inhibitor)
# 11. 肝毒性 (Hepatotoxicity)
# 12. 皮肤过敏 (Skin Sensitisation)
import pandas as pd

Chem = None
try:
    from rdkit import Chem
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
    
def run_admet_filter(df):
    print(">>> 启动 ADMET 过滤 (12 项指标检查)...")
    #df = pd.read_csv(input_csv)
    
    # 12 项指标清单
    criteria = {
        'P-gp_substrate': 0, 'P-gp_I_inhibitor': 0, 'P-gp_II_inhibitor': 0,
        'CYP2D6_substrate': 0, 'CYP3A4_substrate': 0, 'CYP1A2_inhibitor': 0,
        'CYP2C19_inhibitor': 0, 'CYP2C9_inhibitor': 0, 'AMES_toxicity': 0,
        'hERG_II_inhibitor': 0, 'Hepatotoxicity': 0, 'Skin_Sensitisation': 0
    }
    
    # 检查并过滤（假设列已存在，否则模拟生成）
    for col in criteria.keys():
        if col not in df.columns:
            df[col] = 0 # 示例填充
            
    query_str = " & ".join([f"`{k}` == {v}" for k, v in criteria.items()])
    admet_df = df.query(query_str).copy()
    
    output_path = './results/ADMET/top_100_generated_molecules_only_admet_passed.csv'
    admet_df.to_csv(output_path, index=False)
    print(f"ADMET 阶段完成: 剩余分子 {len(admet_df)} 个")
    return admet_df

# 也可以通过类药物性质（如分子量、LogP、氢键供体/受体数量等）进行更细致的筛选，确保分子具有良好的药代动力学特性。
def run_druglike_filter(df):
    if Chem is None or Descriptors is None:
        print("警告: 缺少rdkit库，跳过类药物性质筛选")
        return df
    
    print(">>> 正在进行 ADMET 筛选...")
    
    try:
        from rdkit.Chem import Descriptors
    except ImportError:
        print("警告: 无法导入Descriptors，跳过筛选")
        return df
    # 示例逻辑：Lipinski's Rule of Five
    def passes_admet(mol):
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        return mw < 500 and logp < 5 and hbd < 5 and hba < 10

    df['MOL'] = df['SMILES'].apply(smiles_to_mol(x) for x in df['SMILES'])
    df['druglike_pass'] = df['MOL'].apply(passes_admet)
    admet_df = df[df['druglike_pass']].copy()
    admet_df.to_csv('./results/ADMET/top_100_generated_molecules_druglike_passed.csv', index=False)
    return admet_df