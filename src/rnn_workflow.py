#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
RNN分子生成模块
RNN Molecule Generation Module

功能: 使用RNN生成新的分子结构
"""

import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path

selfies = None
Chem = None
tf = None

try:
    import selfies as sf
    from rdkit import Chem
    import tensorflow as tf
except ImportError:
    pass


def run_rnn_generation(config=None, demo_mode=False):
    """
    运行RNN分子生成
    
    参数:
        config: 配置字典
        demo_mode: 是否为演示模式
    
    返回:
        生成的分子SMILES列表
    """
    print("=" * 60)
    print("步骤1: RNN分子生成")
    print("=" * 60)
    
    if config is None:
        config = {}
    
    num_molecules = config.get('num_molecules', 500)
    temperature = config.get('temperature', 0.7)
    train_data_path = config.get('train_data', 'data/Anticancer_compounds@1994_from_paper.csv')
    model_path = config.get('model_path', 'models/selfies_generator_rnn.keras')
    output_path = config.get('output_path', 'results/generated_smiles.csv')
    
    if demo_mode:
        print("演示模式: 使用预生成的示例分子")
        demo_molecules = [
            "CCOc1ccc(CC(=O)O)cc1",
            "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
            "COc1ccc(CC(=O)O)cc1",
            "CC(C)Cc1ccc(CC(=O)O)cc1",
            "COc1ccc(C(C)C(=O)O)cc1",
            "c1ccc(C(=O)O)cc1",
            "CC(=O)Oc1ccccc1C(=O)O",
            "OC(=O)c1ccccc1O",
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            "OC(=O)Cc1ccccc1Nc2c(Cl)cccc2Cl"
        ]
        
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        pd.DataFrame({'SMILES': demo_molecules}).to_csv(output_path, index=False)
        print(f"演示分子已保存至: {output_path}")
        print(f"生成了 {len(demo_molecules)} 个示例分子")
        return demo_molecules
    
    if tf is None or Chem is None:
        print("错误: 缺少必要的依赖库 (tensorflow, rdkit)")
        print("请安装: pip install tensorflow rdkit selfies")
        return []
    
    try:
        if os.path.exists(model_path):
            print(f"加载预训练模型: {model_path}")
            model = tf.keras.models.load_model(model_path)
            print("模型加载成功")
            
            alphabet = _load_alphabet()
            generated_smiles = _generate_molecules(model, alphabet, num_molecules, temperature)
            
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            pd.DataFrame({'SMILES': generated_smiles}).to_csv(output_path, index=False)
            print(f"生成了 {len(generated_smiles)} 个分子")
            print(f"结果已保存至: {output_path}")
            
            return generated_smiles
        else:
            print(f"警告: 模型文件不存在: {model_path}")
            print("请先训练RNN模型或使用演示模式")
            return []
            
    except Exception as e:
        print(f"RNN生成失败: {str(e)}")
        import traceback
        traceback.print_exc()
        return []


def _load_alphabet():
    """加载SELFIES字母表"""
    default_alphabet = [
        '[C]', '[c]', '[N]', '[n]', '[O]', '[o]', '[S]', '[s]',
        '[F]', '[Cl]', '[Br]', '[I]', '[P]', '[B]',
        '[#C]', '[=C]', '[=N]', '[=O]', '[=S]',
        '[Ring1]', '[Ring2]', '[Ring3]', '[Branch1]', '[Branch2]',
        '[epsilon]'
    ]
    return default_alphabet


def _generate_molecules(model, alphabet, num_molecules, temperature):
    """使用模型生成分子"""
    generated_smiles = []
    token2idx = {tok: idx + 1 for idx, tok in enumerate(alphabet)}
    idx2token = {idx: tok for tok, idx in token2idx.items()}
    
    max_len = 50
    
    for _ in range(num_molecules):
        try:
            seq = [token2idx.get('[C]', 1)]
            
            for _ in range(max_len - 1):
                x = np.array([seq])
                preds = model.predict(x, verbose=0)
                next_token_idx = _sample_with_temperature(preds[0, len(seq)-1], temperature)
                
                if next_token_idx == 0:
                    break
                    
                seq.append(next_token_idx)
            
            selfies_str = ''.join([idx2token.get(idx, '') for idx in seq if idx > 0])
            smiles = sf.decoder(selfies_str)
            
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                canonical = Chem.MolToSmiles(mol)
                if canonical not in generated_smiles:
                    generated_smiles.append(canonical)
                    
        except Exception:
            continue
    
    return generated_smiles


def _sample_with_temperature(preds, temperature):
    """温度采样"""
    preds = np.asarray(preds).astype('float64')
    preds = np.log(preds + 1e-10) / temperature
    exp_preds = np.exp(preds)
    preds = exp_preds / np.sum(exp_preds)
    probas = np.random.multinomial(1, preds, 1)
    return np.argmax(probas)


if __name__ == "__main__":
    smiles_list = run_rnn_generation(demo_mode=True)
    print(f"生成了 {len(smiles_list)} 个分子")
