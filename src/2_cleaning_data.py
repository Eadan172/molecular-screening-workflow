import os,sys
import argparse
import pandas as pd
import multiprocessing as mp
from tqdm import tqdm
import shutil
import zipfile
from rdkit import Chem
from rdkit.Chem import Descriptors

def calculate_mw(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Descriptors.MolWt(mol)
    except:
        pass
    return 0

def data_cleaning_subprocess_for_zipfile(file_path):
    try:
        # 1. 读取数据（增加 low_memory 优化）
        df = pd.read_csv(file_path, encoding='utf-8-sig', low_memory=False)
        
        # 2. 严格单位过滤 (解决 M ml-1 问题)
        pattern = r'^(nM|uM|mM)$'
        mask = df['Standard Units'].str.fullmatch(pattern, case=False, na=False)
        df = df[mask].copy()
        
        # 3. 数值清洗：确保 IC50 是数字且大于 0
        df['Standard Value'] = pd.to_numeric(df['Standard Value'], errors='coerce')
        df = df[df['Standard Value'].notnull() & (df['Standard Value'] > 0)]
        
        # 4. 列筛选与重命名
        cols = ['Molecule ChEMBL ID', 'Smiles','Molecular Weight','Standard Value', 
                'Standard Units','Target ChEMBL ID','Target Name']
        df = df[cols].rename(columns={
            'Molecule ChEMBL ID' : 'ChEMBL ID',
            'Smiles': 'SMILES',
            'Molecular Weight': 'MW',
            'Standard Value' : 'IC50'
        })
        
        # 5. 单位统一化 (可选：全部转为 nM 方便对比)
        # uM * 1000, mM * 1000000
        def convert_to_nm(row):
            unit = str(row['Standard Units']).lower()
            val = row['IC50']
            if unit == 'um': return val * 1000
            if unit == 'mm': return val * 1000000
            return val
        
        df['IC50_nM'] = df.apply(convert_to_nm, axis=1)
        
        return df
    
    except Exception as e:
        print(f"处理文件 {file_path} 时出错: {e}")
        return pd.DataFrame()

def data_cleaning_subprocess_for_api(file_path):
    try:
        # 1. 读取数据（增加 low_memory 优化）
        df = pd.read_csv(file_path, encoding='utf-8-sig', low_memory=False)
        
        # 2. 严格单位过滤 (解决 M ml-1 问题)
        pattern = r'^(nM|uM|mM)$'
        mask = df['standard_units'].str.fullmatch(pattern, case=False, na=False)
        df = df[mask].copy()
        
        # 3. 数值清洗：确保 IC50 是数字且大于 0
        df['standard_value'] = pd.to_numeric(df['standard_value'], errors='coerce')
        df = df[df['standard_value'].notnull() & (df['standard_value'] > 0)]
        
        # 4. 列筛选与重命名
        cols = ['molecule_chembl_id', 'canonical_smiles','standard_value', 
                'standard_units','target_chembl_id','target_pref_name']
        df = df[cols].rename(columns={
            'molecule_chembl_id' : 'ChEMBL ID',
            'canonical_smiles': 'SMILES',
            'standard_value' : 'IC50',
            'target_pref_name' : 'Target Name',
            'standard_units' : 'Standard Units',
            'target_chembl_id' : 'Target ChEMBL ID'
        })
        print(f"处理后的列名: {df.columns.tolist()}")
        # 5.生成分子量列（API数据中没有，需要计算）
        df['MW'] = df['SMILES'].apply(calculate_mw)
        print(len(df['MW'])==len(df))

        # 6. 单位统一化 (可选：全部转为 nM 方便对比)
        # uM * 1000, mM * 1000000
        def convert_to_nm(row):
            unit = str(row['Standard Units']).lower()
            val = row['IC50']
            if unit == 'um': return val * 1000
            if unit == 'mm': return val * 1000000
            return val
        
        df['IC50_nM'] = df.apply(convert_to_nm, axis=1)
        
        # 重新排序列
        df = df[['ChEMBL ID', 'SMILES', 'MW', 'IC50', 'Standard Units', 'Target ChEMBL ID', 'Target Name', 'IC50_nM']]
        return df
    
    except Exception as e:
        print(f"处理文件 {file_path} 时出错: {e}")
        return pd.DataFrame()

# 复制parallel_process_files函数
def parallel_process_files(file_list, final_output, processing_function):
    tmp_dir = './tmp_processing_test'
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # 准备任务参数：只传递文件路径
    tasks = file_list
    print(f"任务列表: {tasks}")

    num_cores = max(1, 3) 
    print(f"启动多核处理，核心数: {num_cores}")

    with mp.Pool(processes=num_cores) as pool:
        results = list(tqdm(pool.imap_unordered(processing_function, tasks), 
                  total=len(tasks), desc="清洗进度"))

    # --- 合并阶段 ---
    print("正在合并数据...")
    print(f"处理结果数量: {len(results)}")
    print(f"第一个结果类型: {type(results[0])}")
    
    # 过滤空的 DataFrame
    non_empty_results = [df for df in results if not df.empty]
    print(f"非空结果数量: {len(non_empty_results)}")
    
    if not non_empty_results:
        print("没有产生任何有效数据。")
        return

    # 合并所有结果
    final_df = pd.concat(non_empty_results, ignore_index=True)
    print(f"合并后的数据形状: {final_df.shape}")
    final_df.to_csv(final_output, index=False, encoding='utf-8-sig')

    
    file_size_gb = os.path.getsize(final_output) / (1024**3)
    if file_size_gb > 1.0:
        zip_name = final_output.replace('.csv', '.zip')
        print(f"文件大小为 {file_size_gb:.2f}GB (超过1GB)，正在压缩为 {zip_name}...")
        with zipfile.ZipFile(zip_name, 'w', zipfile.ZIP_DEFLATED) as zf:
            zf.write(final_output, arcname=os.path.basename(final_output))
        os.remove(final_output) # 删除原大文件
        print("压缩完成，已删除原 CSV。")

    
    shutil.rmtree(tmp_dir)
    print(f"任务结束！清理了临时目录。")

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description='并行处理CHEMBL数据清洗')
    argparser.add_argument('processing_function', choices=['zip', 'api'], help='选择处理函数 zip:处理直接从CHEMBL下载的压缩包 api:处理通过API获取的数据\n')
    args = argparser.parse_args()

    
    batch_dir = './chembl_split_data_zip'
    
    if not os.path.exists(batch_dir):
        print(f"目录 {batch_dir} 不存在")
    else:
        all_batches = [os.path.join(batch_dir, f) for f in os.listdir(batch_dir) if f.endswith('.csv')]
        print(f"找到 {len(all_batches)} 个CSV文件")
        
        
        all_batches.sort()
        print(f"文件列表: {all_batches}")
    
    
    if sys.argv[1] == 'zip':
        parallel_process_files(all_batches, './ChEMBL_Cleaned_Finall_test_zipfile.csv', data_cleaning_subprocess_for_zipfile)
    elif sys.argv[1] == 'api':
        parallel_process_files(all_batches, './ChEMBL_Cleaned_Finall_test_api.csv', data_cleaning_subprocess_for_api)
