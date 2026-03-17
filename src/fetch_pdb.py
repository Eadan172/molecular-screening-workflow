import os
import json
import requests
from Bio import PDB

def fetch_pdb_ids(target_name, organism, max_resolution):
    """从 RCSB 搜索符合条件的 PDB ID 列表"""
    url = "https://search.rcsb.org/rcsbsearch/v2/query?json="
    
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                        "operator": "exact_match",
                        "value": target_name
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name",
                        "operator": "exact_match",
                        "value": organism
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "less_or_equal",
                        "value": max_resolution
                    }
                }
            ]
        },
        "return_type": "entry"
    }
    
    try:
        response = requests.get(url + json.dumps(query))
        response.raise_for_status() # 检查 HTTP 响应状态
        response_data = response.json()
        
        if 'result_set' in response_data:
            pdb_ids = [entry['identifier'] for entry in response_data['result_set']]
            return pdb_ids
        else:
            print(f"提示: 未找到符合条件的结构。")
            return []
    except Exception as e:
        print(f"搜索出错: {e}")
        return []

def download_and_rename_pdb(pdb_ids, target_name):
    """下载 PDB 文件并将 .ent 格式统一修改为 .pdb 格式"""
    if not pdb_ids:
        print("没有可下载的 ID。")
        return

    # 创建保存目录
    save_dir = f"./pdb_files_{target_name}"
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
        print(f"已创建目录: {save_dir}")

    pdbl = PDB.PDBList()
    
    for pdb_id in pdb_ids:
        # 1. 执行下载 (默认会下载到 save_dir 下的 pdb 文件夹中，并命名为 pdbXXXX.ent)
        # 注意: retrieve_pdb_file 返回的是下载后的文件完整路径
        downloaded_path = pdbl.retrieve_pdb_file(
            pdb_id, 
            file_format="pdb", 
            pdir=save_dir, 
            overwrite=True
        )
        
        # 2. 检查文件是否存在并重命名
        if os.path.exists(downloaded_path):
            # 构造新的文件名，例如: 1abc.pdb
            new_filename = f"{pdb_id.lower()}.pdb"
            new_path = os.path.join(save_dir, new_filename)
            
            # 执行重命名操作 (从 .../pdbXXXX.ent 到 .../XXXX.pdb)
            os.rename(downloaded_path, new_path)
            print(f"成功转换: {new_filename}")
        else:
            print(f"警告: 无法找到已下载的文件 {pdb_id}")

# --- 主程序运行 ---

# 1. 查询参数
target_gene = "AKT1" 
target_organism = "Homo sapiens"
resolution_cutoff = 2.0

# 2. 执行搜索
print(f"正在搜索 {target_gene} ({target_organism})，分辨率 <= {resolution_cutoff}Å...")
found_ids = fetch_pdb_ids(target_gene, target_organism, resolution_cutoff)
print(f"找到 {len(found_ids)} 个符合条件的 PDB IDs: {found_ids}")

# 3. 执行下载与重命名
download_and_rename_pdb(found_ids, target_gene)

print("\n任务完成！文件已存放在相应的文件夹中。")