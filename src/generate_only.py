import pandas as pd
import numpy as np
import tensorflow as tf
import selfies as sf
from rdkit import Chem
import random
import os

print("Loading data...")
df1994 = pd.read_csv('./csv_input/Anticancer_compounds@1994_from_paper.csv', encoding='latin1')
print('Data loaded:', str(df1994.shape))

def valid_smiles(s):
    try:
        return Chem.MolToSmiles(Chem.MolFromSmiles(s))
    except:
        return None

df1994['SMILES'] = df1994['SMILES'].apply(valid_smiles)
df1994 = df1994.dropna(subset=['SMILES']).reset_index(drop=True)
print('Valid molecules:', len(df1994))

df1994['SELFIES'] = [sf.encoder(s) for s in df1994['SMILES']]
df1994 = df1994.dropna(subset=['SELFIES'])
df1994 = df1994[df1994['SELFIES'].apply(lambda x: isinstance(x, str) and len(x) > 0)]
df1994 = df1994.reset_index(drop=True)
print('SELFIES cleaning complete.')

alphabet = list(sf.get_alphabet_from_selfies(df1994['SELFIES']))
alphabet.append('.')
token2idx = {tok: idx + 1 for idx, tok in enumerate(alphabet)}
idx2token = {idx: tok for tok, idx in token2idx.items()}
max_len = max(len(list(sf.split_selfies(s))) for s in df1994['SELFIES'])
print("Max SELFIES length:", max_len)
print("Alphabet size:", len(alphabet))

print("Loading model...")
if os.path.exists('./models/selfies_generator_rnn.keras'):
    model = tf.keras.models.load_model('./models/selfies_generator_rnn.keras')
    print("Model loaded successfully from selfies_generator_rnn.keras")
else:
    print("Model file not found!")
    exit(1)

def sample_selfies(n, temperature=0.7, top_k=50, model=None):
    new_selfies = []
    padding_token_idx = token2idx.get('.', 0)
    
    print(f"Generating {n} molecules...")

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
            
            probs_for_next_token = np.log(probs_for_next_token + 1e-8) / temperature
            probs_for_next_token = np.exp(probs_for_next_token) / np.sum(np.exp(probs_for_next_token))
            
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
                        if 10 <= len(smi) <= 100:
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

print("Generating molecules...")
generated_smiles = sample_selfies(500, temperature=0.7, top_k=50, model=model)

all_generated_df = pd.DataFrame(generated_smiles, columns=['SMILES'])
all_generated_df.to_csv('./results/all_generated_molecules.csv', index=False)
print(f"All generated molecules saved to ./results/all_generated_molecules.csv")

generated_smiles = generated_smiles[:100]
smiles_df = pd.DataFrame(generated_smiles, columns=['SMILES'])
output_filename = './results/generated_smiles_from_1994mols.csv'
smiles_df.to_csv(output_filename, index=False)
print(f'Top 100 molecules saved to {output_filename}')
print('Done!')
