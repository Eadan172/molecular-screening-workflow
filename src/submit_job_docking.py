#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
分子对接模块
Molecular Docking Module

功能: 提交分子对接任务
"""

import os
import pandas as pd
import json
import requests

Chem = None
AllChem = None
SDWriter = None
rdmolfiles = None

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, SDWriter
    from rdkit.Chem import rdmolfiles
except ImportError:
    pass

try:
    from requests_toolbelt.multipart.encoder import MultipartEncoder
except ImportError:
    MultipartEncoder = None


def csv_to_sdf(csv_file, sdf_file):
    """
    将 CSV 转换为 SDF 文件。

    参数:
        csv_file (str): 输入 CSV 文件路径。
        sdf_file (str): 输出 SDF 文件路径。
    """
    if Chem is None or SDWriter is None:
        print("警告: 缺少rdkit库，无法转换为SDF")
        return False
    
    df = pd.read_csv(csv_file)
    writer = SDWriter(sdf_file)

    for index, row in df.iterrows():
        smiles = row['SMILES']
        mol = Chem.MolFromSmiles(smiles)
        
        if mol:
            AllChem.Compute2DCoords(mol)
            for col in df.columns:
                mol.SetProp(col, str(row[col]))
            writer.write(mol)
        else:
            print(f"无法处理 SMILES: {smiles} (行索引: {index})")

    writer.close()
    print(f"SDF 文件'{sdf_file}'已生成！")
    return True


def csv_to_mol2(csv_file, output_folder):
    """
    将 CSV 文件中的 SMILES 转换为 MOL2 文件。
    
    参数:
        csv_file (str): 输入 CSV 文件路径。
        output_folder (str): 存储 MOL2 文件的目标文件夹路径。
    """
    if Chem is None or rdmolfiles is None:
        print("警告: 缺少rdkit库，无法转换为MOL2")
        return False
    
    df = pd.read_csv(csv_file)
    os.makedirs(output_folder, exist_ok=True)

    for index, row in df.iterrows():
        smiles = row["SMILES"]
        name = row["Name"] if "Name" in row else f"Molecule_{index}"
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"无法解析 SMILES: {smiles} (行索引: {index})")
            continue
        
        mol_with_h = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_with_h, randomSeed=42)
        AllChem.UFFOptimizeMolecule(mol_with_h)

        mol2_file = f"{output_folder}/{name}.mol2"
        
        try:
            rdmolfiles.MolToMol2File(mol_with_h, mol2_file)
            print(f"MOL2 文件已保存: {mol2_file}")
        except Exception as e:
            print(f"无法保存 MOL2 文件: {mol2_file}, 错误: {e}")
    
    return True


def submit_docking_job(docking_engine=['vina'], demo_mode=False):
    """
    提交分子对接任务
    
    参数:
        docking_engine: 对接引擎列表
        demo_mode: 是否为演示模式
    
    返回:
        任务ID或结果
    """
    print("=" * 60)
    print("步骤4: 分子对接")
    print("=" * 60)
    
    if demo_mode:
        print("演示模式: 跳过实际对接")
        demo_result = pd.DataFrame({
            'SMILES': ['CCOc1ccc(CC(=O)O)cc1'],
            'binding_affinity': [-8.5],
            'rmsd_lb': [0.0],
            'rmsd_ub': [1.2]
        })
        os.makedirs('./results/docking', exist_ok=True)
        demo_result.to_csv('./results/docking/docking_result.csv', index=False)
        print("演示对接结果已保存")
        return demo_result
    
    if MultipartEncoder is None:
        print("警告: 缺少requests_toolbelt库，无法提交在线对接任务")
        print("请安装: pip install requests-toolbelt")
        return None
    
    docker_dict = {
        "gnina": "https://neurosnap.ai/api/job/submit/GNINA?note=my job description",
        "diffdock-l": "https://neurosnap.ai/api/job/submit/DiffDock-L",
        "autodock-vina": "https://neurosnap.ai/api/job/submit/AutoDock Vina (smina)?note=my job description"
    }
    
    target_name = 'AKT1'
    pdb_dir = './pdb_files_' + target_name
    
    if not os.path.exists(pdb_dir):
        print(f"警告: PDB目录不存在: {pdb_dir}")
        return None
    
    pdb_files = os.listdir(pdb_dir)
    if not pdb_files:
        print("警告: PDB目录为空")
        return None
    
    sdf_file = "./results/ADMET/all_generated_molecules_only_admet_passed.sdf"
    if not os.path.exists(sdf_file):
        print(f"警告: 配体SDF文件不存在: {sdf_file}")
        return None
    
    try:
        multipart_data = MultipartEncoder(
            fields={
                "Input Receptor": ("structure.pdb", open(os.path.join(pdb_dir, pdb_files[0]), "rb")),
                "Input Ligand": json.dumps([{"data": open(sdf_file).read(), "type": "sdf"}]),
            }
        )
        
        YOUR_API_KEY = os.environ.get("NEUROSNAP_API_KEY", "your_api_key_here")
        
        engine = docking_engine[0] if docking_engine else 'autodock-vina'
        url = docker_dict.get(engine, docker_dict['autodock-vina'])
        
        r = requests.post(
            url,
            headers={
                "X-API-KEY": YOUR_API_KEY,
                "Content-Type": multipart_data.content_type,
            },
            data=multipart_data,
        )
        
        new_job_id = r.json()
        print(f"对接任务已提交，任务ID: {new_job_id}")
        return new_job_id
        
    except Exception as e:
        print(f"提交对接任务失败: {e}")
        return None


if __name__ == "__main__":
    result = submit_docking_job(demo_mode=True)
    print(result)
