#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
SA Score评估模块
Synthetic Accessibility Score Evaluation Module

功能: 计算分子的合成可及性评分
"""

import os
import sys
import pandas as pd
import numpy as np

Chem = None
sascorer = None

try:
    from rdkit import Chem
    from rdkit.Chem import RDConfig
    sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
    import sascorer
except ImportError:
    pass


def calculate_sa_score(smiles):
    """
    计算单个分子的SA Score
    
    参数:
        smiles: 分子SMILES字符串
    
    返回:
        SA Score (1-10)
    """
    if sascorer is None:
        return None
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return sascorer.calculateScore(mol)
    except Exception:
        return None


def run_sa_filter(df, threshold=6.0, output_path='results/SA/sa_filtered.csv'):
    """
    运行SA Score筛选
    
    参数:
        df: 包含SMILES列的数据框
        threshold: SA Score阈值
        output_path: 输出文件路径
    
    返回:
        筛选后的数据框
    """
    print("=" * 60)
    print("步骤5: SA Score评估")
    print("=" * 60)
    
    if sascorer is None:
        print("警告: SA Score模块不可用，跳过此步骤")
        return df
    
    df = df.copy()
    
    if 'SMILES' not in df.columns:
        print("错误: 数据框缺少SMILES列")
        return df
    
    print(f"计算 {len(df)} 个分子的SA Score...")
    
    df['mol'] = df['SMILES'].apply(lambda x: Chem.MolFromSmiles(x) if pd.notna(x) else None)
    df['sa_score'] = df['mol'].apply(lambda x: sascorer.calculateScore(x) if x is not None else None)
    
    df = df.dropna(subset=['sa_score'])
    
    print(f"SA Score范围: {df['sa_score'].min():.2f} - {df['sa_score'].max():.2f}")
    print(f"SA Score均值: {df['sa_score'].mean():.2f}")
    
    filtered_df = df[df['sa_score'] <= threshold].copy()
    filtered_df = filtered_df.sort_values('sa_score')
    
    if 'mol' in filtered_df.columns:
        filtered_df = filtered_df.drop(columns=['mol'])
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    filtered_df.to_csv(output_path, index=False)
    
    print(f"通过SA Score筛选的分子: {len(filtered_df)} / {len(df)}")
    print(f"结果已保存至: {output_path}")
    
    return filtered_df


def get_sa_score_statistics(df):
    """
    获取SA Score统计信息
    
    参数:
        df: 包含sa_score列的数据框
    
    返回:
        统计字典
    """
    if 'sa_score' not in df.columns:
        return None
    
    return {
        'min': df['sa_score'].min(),
        'max': df['sa_score'].max(),
        'mean': df['sa_score'].mean(),
        'median': df['sa_score'].median(),
        'std': df['sa_score'].std()
    }


if __name__ == "__main__":
    test_smiles = "CCOc1ccc(CC(=O)O)cc1"
    score = calculate_sa_score(test_smiles)
    print(f"测试分子 SA Score: {score}")
