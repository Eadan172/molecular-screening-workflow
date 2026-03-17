import pandas as pd
# 确保CSV Files文件夹在代码所在目录下
df1994 = pd.read_csv('./csv_input/Anticancer_compounds@1994_from_paper.csv', encoding='latin1')
print('Data loaded:', str(df1994.shape))

import rdkit
from rdkit import Chem
def valid_smiles(s):
    try:
        return Chem.MolToSmiles(Chem.MolFromSmiles(s))
    except:
        return None

df1994['SMILES'] = df1994['SMILES'].apply(valid_smiles)
df1994 = df1994.dropna(subset=['SMILES']).reset_index(drop=True)
print('Valid molecules:', len(df1994))

import selfies as sf
df1994['SELFIES'] = [sf.encoder(s) for s in df1994['SMILES']]

df1994 = df1994.dropna(subset=['SELFIES'])
df1994 = df1994[df1994['SELFIES'].apply(lambda x: isinstance(x, str) and len(x) > 0)]
df1994 = df1994.reset_index(drop=True)
print('SELFIES cleaning complete.')

alphabet = list(sf.get_alphabet_from_selfies(df1994['SELFIES']))
print(alphabet)
print(len(alphabet))
alphabet.append('.')

token2idx = {tok: idx + 1 for idx, tok in enumerate(alphabet)}
idx2token = {idx: tok for tok, idx in token2idx.items()} 

print('token2idx',token2idx)
print('idx2token',idx2token)

max_len = max(len(list(sf.split_selfies(s))) for s in df1994['SELFIES'])
print("Max SELFIES length:", max_len)
print("Alphabet size:", len(alphabet))

import numpy as np
# 方法1：selfies字符串转为索引列表,并进行padding
def encode_sequence(s):
    toks = sf.split_selfies(s)
    seq = [token2idx[tok] for tok in toks]
    # Pad to max_len
    return seq + [0] * (max_len - len(seq))

seqs = np.array([encode_sequence(s) for s in df1994['SELFIES']], dtype=np.int32)

# 方法2：使用Keras的pad_sequences（需要先安装TensorFlow）
from tensorflow.keras.preprocessing.sequence import pad_sequences as tf_pad

def encode_sequence_tf(s, token2idx):
    """将字符串序列转换为数字ID序列"""
    # 逐个字符映射，未知字符用<UNK>的ID,这里是0
    toks = sf.split_selfies(s)
    return [token2idx.get(tok, 0) for tok in toks]

# 编码所有序列
encoded_sequences = [encode_sequence_tf(s, token2idx) for s in df1994['SELFIES']]

# 直接调用内置函数
tf_padded = tf_pad(encoded_sequences, maxlen=max_len, padding="post", value=0)

# 方法3：使用PyTorch的pad_sequence（需要先转换为tensor）
import torch
from torch.nn.utils.rnn import pad_sequence

def encode_sequence_tf(s, token2idx):
    """将字符串序列转换为数字ID序列"""
    # 逐个字符映射，未知字符用<UNK>的ID
    toks = sf.split_selfies(s)
    return [token2idx.get(tok, 0) for tok in toks]

# 编码所有序列
encoded_sequences = [encode_sequence_tf(s, token2idx) for s in df1994['SELFIES']]

# 将编码序列转为tensor
tensor_seqs = [torch.tensor(seq) for seq in encoded_sequences]
# 填充（默认填充到最长序列，padding_value是填充符ID）
torch_padded = pad_sequence(tensor_seqs, batch_first=True, padding_value=0) # 用0作为padding值
print("PyTorch填充结果：\n", torch_padded)

X_gen = seqs[:, :-1] #or X_gen = tf_padded[:, :-1] or X_gen = torch_padded[:, :-1]
y_gen = seqs[:, 1:] #or y_gen = tf_padded[:, 1:] or y_gen = torch_padded[:, 1:]
print(seqs.shape,X_gen.shape,y_gen.shape)

# 定义RNN模型
import tensorflow as tf
vocab_size = len(token2idx) + 1

input_layer = tf.keras.Input(shape=(max_len - 1,))

embedding_layer = tf.keras.layers.Embedding(vocab_size, 128, mask_zero=True)(input_layer)
dropout1 = tf.keras.layers.Dropout(0.2)(embedding_layer)
lstm1 = tf.keras.layers.LSTM(256, return_sequences=True, recurrent_dropout=0.2)(dropout1)
dropout2 = tf.keras.layers.Dropout(0.2)(lstm1)
lstm2 = tf.keras.layers.LSTM(256, return_sequences=True, recurrent_dropout=0.2)(dropout2)
dropout3 = tf.keras.layers.Dropout(0.2)(lstm2)
lstm3 = tf.keras.layers.LSTM(128, return_sequences=True, recurrent_dropout=0.2)(dropout3)

output_layer = tf.keras.layers.Dense(vocab_size, activation='softmax')(lstm3)

model = tf.keras.Model(inputs=input_layer, outputs=output_layer)

optimizer = tf.keras.optimizers.Adam(learning_rate=0.001)
model.compile(optimizer=optimizer, loss='categorical_crossentropy',metrics=['accuracy','precision'])

reduce_lr = tf.keras.callbacks.ReduceLROnPlateau(
    monitor='loss', 
    factor=0.5, 
    patience=3, 
    verbose=1
)

early_stopping = tf.keras.callbacks.EarlyStopping(
    monitor='loss',
    patience=10,
    verbose=1,
    restore_best_weights=True
)


from tensorflow.keras.utils import to_categorical
def data_generator(X, y, batch_size=64, num_classes=vocab_size):
    while True:
        for i in range(0, len(X), batch_size):
            X_batch = X[i:i+batch_size]
            y_batch = y[i:i+batch_size]
            y_batch = to_categorical(y_batch, num_classes=num_classes)
            
            yield X_batch, y_batch

train_gen = data_generator(X_gen, y_gen, batch_size=32)

steps_per_epoch = max(1, len(X_gen) // 32) # 确保至少有1步,//的意思是整数除法，确保得到一个整数结果，避免出现小数步数导致训练出错。

# 测试模型输出维度是否正确
dummy_input = np.zeros((1, max_len - 1), dtype=np.int32) 

print("✅ Dummy input shape:", dummy_input.shape) 
print("✅ Model output shape:", model(dummy_input).shape) 


history = model.fit(
    train_gen,
    steps_per_epoch=steps_per_epoch,
    epochs=100,
    verbose=1,
    callbacks=[reduce_lr, early_stopping]
)

import json
from datetime import datetime

rnn_performance = {
    'model_type': 'LSTM (3-layer)',
    'architecture': {
        'lstm_layers': [256, 256, 128],
        'dropout': 0.2,
        'vocab_size': vocab_size,
        'max_sequence_length': max_len
    },
    'training': {
        'epochs_trained': len(history.history['loss']),
        'final_loss': round(float(history.history['loss'][-1]), 4),
        'final_accuracy': round(float(history.history['accuracy'][-1]), 4),
        'final_precision': round(float(history.history['precision'][-1]), 4) if 'precision' in history.history else None,
        'optimizer': 'Adam',
        'learning_rate': 0.001
    },
    'data': {
        'training_samples': len(X_gen),
        'alphabet_size': len(alphabet)
    },
    'training_date': datetime.now().isoformat()
}

os.makedirs('./models', exist_ok=True)
with open('./models/rnn_performance.json', 'w', encoding='utf-8') as f:
    json.dump(rnn_performance, f, indent=2, ensure_ascii=False)
print('模型性能记录已保存: ./models/rnn_performance.json')

import random


def sample_selfies(n, temperature=0.7, top_k=50, model=None):
    new_selfies = []
    
    padding_token_idx = token2idx.get('.', 0) 

    for _ in range(n): 
        
        valid_tokens = [idx for token, idx in token2idx.items() if token != '.'] 
        
        if not valid_tokens:
            continue

        seq = [random.choice(valid_tokens)]
        
        
        for i in range(max_len - 2):
            current_seq = seq[:max_len - 1]
            
            current_padded_seq = current_seq + [padding_token_idx] * (max_len - 1 - len(current_seq))
            
            current_padded_seq = np.array([current_padded_seq], dtype=np.int32)
            
            probs = model.predict(current_padded_seq, verbose=0)[0]
            
            probs_for_next_token = probs[len(current_seq) - 1]
            
            # 应用温度参数
            probs_for_next_token = np.log(probs_for_next_token + 1e-8) / temperature
            probs_for_next_token = np.exp(probs_for_next_token) / np.sum(np.exp(probs_for_next_token))
            
            # 应用top-k采样
            top_k_indices = np.argsort(probs_for_next_token)[-top_k:]
            top_k_probs = probs_for_next_token[top_k_indices]
            top_k_probs = top_k_probs / np.sum(top_k_probs)
            
            next_tok_idx = np.random.choice(top_k_indices, p=top_k_probs)
            
            if idx2token.get(next_tok_idx, '') == '.' or next_tok_idx == padding_token_idx:
                 break
            

            seq.append(next_tok_idx)
            

        toks = [idx2token.get(i, '') for i in seq if idx2token.get(i, '') != '.']
        
        generated_selfie = ''.join(toks)
        
        if generated_selfie: 
            try:
                
                smi = sf.decoder(generated_selfie)
                if smi is not None:
                     mol = Chem.MolFromSmiles(smi)
                     if mol is not None:
                        # 检查分子大小和有效性
                        if 10 <= len(smi) <= 100:  # 合理的分子大小
                            new_selfies.append(generated_selfie)

            except Exception as e:
                pass
      
    generated_smiles = []
    for s in new_selfies:
        try:
            smi = sf.decoder(s)
            if smi is not None:
                mol = Chem.MolFromSmiles(smi)
                if mol is not None:
                    generated_smiles.append(smi)
        except:
            pass


    print('Valid generated SMILES:', len(generated_smiles))
    return generated_smiles

import tensorflow as tf

# 保存模型为正确的格式
model.save('./models/selfies_generator_rnn.keras')
print('模型已保存为 .keras 格式')

# 加载训练好的模型
loaded_model = tf.keras.models.load_model('./models/selfies_generator_rnn.keras')
print("Model loaded successfully.")
# 生成更多分子以确保得到100个有效分子
generated_smiles = sample_selfies(500, temperature=0.7, top_k=50, model=loaded_model)
# 确保只保留前100个分子
generated_smiles = generated_smiles[:100]

smiles_df = pd.DataFrame(generated_smiles, columns=['SMILES'])
output_filename = './results/generated_smiles_from_1994mols.csv'
smiles_df.to_csv(output_filename, index=False)

print(f'Generated SMILES saved to {output_filename} and available for download.')