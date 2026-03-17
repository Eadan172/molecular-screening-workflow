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
<<<<<<< HEAD
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
=======
            try:
                model = tf.keras.models.load_model(model_path)
                print("模型加载成功")
            except Exception as load_error:
                print(f"模型加载失败: {load_error}")
                print("模型文件格式不兼容，将重新训练...")
                model = None
            
            if model is not None:
                alphabet = _load_alphabet()
                generated_smiles = _generate_molecules(model, alphabet, num_molecules, temperature)
                
                os.makedirs(os.path.dirname(output_path), exist_ok=True)
                pd.DataFrame({'SMILES': generated_smiles}).to_csv(output_path, index=False)
                print(f"生成了 {len(generated_smiles)} 个分子")
                print(f"结果已保存至: {output_path}")
                
                return generated_smiles
        
        print("开始训练新的RNN模型...")
        return _train_and_generate(train_data_path, model_path, output_path, num_molecules, temperature)
>>>>>>> 518afaa (2rd commit: modify RNN QSAR model)
            
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


<<<<<<< HEAD
=======
def _train_and_generate(train_data_path, model_path, output_path, num_molecules, temperature):
    """训练RNN模型并生成分子 - 使用优化架构"""
    import json
    from datetime import datetime
    
    print(f"加载训练数据: {train_data_path}")
    df = pd.read_csv(train_data_path, encoding='latin1')
    df = df.dropna(subset=['SMILES']).reset_index(drop=True)
    
    df['SMILES'] = df['SMILES'].apply(lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s)) if Chem.MolFromSmiles(s) else None)
    df = df.dropna(subset=['SMILES']).reset_index(drop=True)
    print(f"有效分子数: {len(df)}")
    
    df['SELFIES'] = [sf.encoder(s) for s in df['SMILES']]
    df = df.dropna(subset=['SELFIES']).reset_index(drop=True)
    
    alphabet = list(sf.get_alphabet_from_selfies(df['SELFIES']))
    alphabet.append('.')
    print(f"字母表大小: {len(alphabet)}")
    
    token2idx = {tok: idx + 1 for idx, tok in enumerate(alphabet)}
    idx2token = {idx: tok for tok, idx in token2idx.items()}
    vocab_size = len(alphabet) + 1
    
    max_len = max(len(list(sf.split_selfies(s))) for s in df['SELFIES']) + 5
    print(f"最大序列长度: {max_len}")
    
    X_gen = []
    y_gen = []
    for selfie in df['SELFIES']:
        tokens = list(sf.split_selfies(selfie))
        token_ids = [token2idx.get(t, 0) for t in tokens]
        for i in range(len(token_ids) - 1):
            X_gen.append(token_ids[:i+1])
            y_gen.append(token_ids[i+1])
    
    from tensorflow.keras.preprocessing.sequence import pad_sequences
    X_gen = pad_sequences(X_gen, maxlen=max_len-1, padding='post')
    y_gen = np.array(y_gen)
    
    print(f"训练样本数: {len(X_gen)}")
    
    # 使用优化的Bi-LSTM架构
    model = tf.keras.Sequential([
        tf.keras.layers.Embedding(vocab_size, 512, input_length=max_len-1),
        tf.keras.layers.Bidirectional(tf.keras.layers.LSTM(512, return_sequences=True, dropout=0.1)),
        tf.keras.layers.LayerNormalization(),
        tf.keras.layers.Bidirectional(tf.keras.layers.LSTM(512, return_sequences=True, dropout=0.1)),
        tf.keras.layers.LayerNormalization(),
        tf.keras.layers.Bidirectional(tf.keras.layers.LSTM(256, return_sequences=False, dropout=0.1)),
        tf.keras.layers.LayerNormalization(),
        tf.keras.layers.Dense(512, activation='relu'),
        tf.keras.layers.Dropout(0.2),
        tf.keras.layers.Dense(vocab_size, activation='softmax')
    ])
    
    model.compile(
        optimizer=tf.keras.optimizers.Adam(learning_rate=0.0005),
        loss='sparse_categorical_crossentropy',
        metrics=['accuracy']
    )
    
    model.summary()
    
    print("开始训练优化模型...")
    history = model.fit(
        X_gen, y_gen,
        epochs=200,
        batch_size=128,
        validation_split=0.1,
        callbacks=[
            tf.keras.callbacks.EarlyStopping(
                monitor='accuracy', 
                patience=10, 
                restore_best_weights=True,
                min_delta=0.001
            ),
            tf.keras.callbacks.ReduceLROnPlateau(
                monitor='accuracy', 
                factor=0.5, 
                patience=5, 
                min_lr=1e-6
            ),
            tf.keras.callbacks.ModelCheckpoint(
                model_path.replace('.keras', '_best.keras'),
                monitor='accuracy',
                save_best_only=True
            )
        ],
        verbose=1
    )
    
    os.makedirs(os.path.dirname(model_path), exist_ok=True)
    model.save(model_path)
    print(f"模型已保存: {model_path}")
    
    final_acc = float(history.history['accuracy'][-1])
    final_val_acc = float(history.history['val_accuracy'][-1]) if 'val_accuracy' in history.history else final_acc
    
    rnn_performance = {
        'model_type': 'Optimized Bi-LSTM (3-layer)',
        'architecture': {
            'lstm_layers': [512, 512, 256],
            'bidirectional': True,
            'layer_normalization': True,
            'embedding_dim': 512,
            'dropout': 0.1,
            'vocab_size': vocab_size,
            'max_sequence_length': max_len
        },
        'training': {
            'epochs_trained': len(history.history['loss']),
            'final_loss': round(float(history.history['loss'][-1]), 4),
            'final_accuracy': round(final_acc, 4),
            'final_val_accuracy': round(final_val_acc, 4),
            'optimizer': 'Adam',
            'initial_learning_rate': 0.0005,
            'batch_size': 128
        },
        'data': {
            'training_samples': len(X_gen),
            'alphabet_size': len(alphabet),
            'unique_molecules': len(df)
        },
        'target_met': final_acc >= 0.95,
        'training_date': datetime.now().isoformat()
    }
    
    perf_path = model_path.replace('.keras', '_performance.json')
    with open(perf_path, 'w', encoding='utf-8') as f:
        json.dump(rnn_performance, f, indent=2, ensure_ascii=False)
    print(f"模型性能记录已保存: {perf_path}")
    
    print("开始生成分子...")
    generated_smiles = _generate_molecules(model, alphabet, num_molecules, temperature)
    
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    pd.DataFrame({'SMILES': generated_smiles}).to_csv(output_path, index=False)
    print(f"生成了 {len(generated_smiles)} 个分子")
    print(f"结果已保存至: {output_path}")
    
    return generated_smiles


>>>>>>> 518afaa (2rd commit: modify RNN QSAR model)
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
