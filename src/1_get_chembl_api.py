#!/usr/bin/env python3
"""
Batch-download ChEMBL molecules using public API that satisfy Lipinski's rule of five (with custom cuts),
returning only IUPAC name, mol weight, IC50 value, SMILES and target name.

Usage:
  python download_chembl_api.py  # Download all targets
  python download_chembl_api.py AKT1 IKKB bcl-2 EGFR p53  # Download specific targets
"""
import requests
import pandas as pd
import time
import os
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

# 配置
BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
PAGE_SIZE = 100       # API 单页上限通常为 100
MAX_WORKERS = 8       # 增加线程数以提高并发度
BATCH_SIZE = 16        # 每次提交的任务批次大小
OUTPUT_DIR = "chembl_split_data_api"
FINAL_ZIP = "chembl_split_data_api.zip"

# 创建会话以保持连接
session = requests.Session()
session.headers.update({
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
})


def get_last_downloaded_page():
    """获取最后下载的页面号"""
    if not os.path.exists(OUTPUT_DIR):
        return 0
    
    max_page = 0
    for file in os.listdir(OUTPUT_DIR):
        if file.startswith('batch_') and file.endswith('.csv'):
            try:
                page_num = int(file.split('_')[1].split('.')[0])
                if page_num > max_page:
                    max_page = page_num
            except:
                pass
    return max_page


if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

def fetch_and_save_batch(page, query_params):
    """下载并直接保存为原始分片文件"""
    params = {
        "offset": (page - 1) * PAGE_SIZE,
        # 计算偏移量，以确保每页数据不重复
        "limit": PAGE_SIZE,
        # API 每页返回的记录数
        **query_params
        # 将基础查询参数合并到请求参数中
    }
    
    for retry in range(3): 
        # 最多重试 3 次以处理临时网络问题，retry 机制可以帮助提高下载的稳定性
        try:
            # 使用会话提高性能
            response = session.get(BASE_URL, params=params, timeout=60, stream=True)
            if response.status_code == 200: 
                # 成功响应
                # 直接解析 JSON，避免中间步骤
                data = response.json().get('activities', [])
                if not data:
                    return None, page # 标记结束
                
                # 直接转为 DataFrame 并保存分片
                df = pd.DataFrame(data)
                batch_file = os.path.join(OUTPUT_DIR, f"batch_{page}.csv")
                # 使用更快的 CSV 写入方式
                df.to_csv(batch_file, index=False, encoding='utf-8-sig', chunksize=1000)
                return len(data), page
            elif response.status_code == 404: 
                # 如果页面不存在，说明已经没有更多数据了
                return None, page
        except Exception as e: 
            # e是捕获到的异常对象，可以包含网络错误、超时等信息
            # 更智能的错误处理
            if 'timeout' in str(e):
                time.sleep(5)  # 超时使用更长的等待时间
            else:
                time.sleep(2 ** retry) 
                # 2的指数级退避，增加等待时间以减少服务器压力
    return 0, page

def main(target_name=None):
    # download_raw_data函数重命名为main，便于直接运行
    # 1. 设置 API 级别的基础过滤（减少无效数据下载）
    query_params = {
        "standard_type": "IC50",
        "target_organism": "Homo sapiens",
    }
    if target_name:
        query_params["target_pref_name__iexact"] = target_name
        # target_pref_name__iexact: 精确匹配目标名称
        # target_pref_name__iexact不存在于chembl直接下载的数据中

    # 获取已下载的最后页面
    last_page = get_last_downloaded_page() #整数类型
    current_page = last_page + 1 if last_page > 0 else 1
    # else 1：从第一页开始下载
    
    if last_page > 0:
        print(f"🔄 检测到已下载至页面 {last_page}，从页面 {current_page} 继续")
    else:
        print(f"🚀 开始并发下载 ChEMBL 原始数据...")
    print(f"配置: 线程数={MAX_WORKERS}, 批次大小={BATCH_SIZE}, 页大小={PAGE_SIZE}")
    
    stop_signal = False
    total_downloaded = 0
    
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:# ThreadPoolExecutor函数
        pbar = tqdm(desc="下载页数")
        
        while not stop_signal:
            # 每次提交一批任务，使用更大的批次大小
            batch_pages = list(range(current_page, current_page + BATCH_SIZE))
            futures = {executor.submit(fetch_and_save_batch, p, query_params): p for p in batch_pages}
            '''
            futures 是一个字典，键是提交的任务（Future对象），
            值是对应的页面号。
            这种结构允许我们在处理完成的任务时，知道它对应的是哪个页面，
            从而更好地管理下载进度和错误处理。
            query_params 是一个字典，这意味着这个字典可能包含了一些键值对，
            用于筛选或定义在 fetch_and_save_batch 函数中的具体操作。
            这种参数通常用于：

            指定查询条件：例如在数据库查询或API请求中，指定需要获取的数据范围或条件。

            配置选项：设定批量操作的选项，例如限制结果数量、排序方式等。

            动态过滤：通过不同的参数组合，动态调整任务的执行逻辑。
            '''
            for future in as_completed(futures):
                count, p_num = future.result()
                if count is None:
                    stop_signal = True
                else:
                    total_downloaded += count
                    pbar.update(1)
            
            current_page += BATCH_SIZE
            # 动态调整延迟，根据任务完成情况
            time.sleep(0.2) # 减少延迟以提高速度

    # 2. 自动打包成 ZIP
    print(f"\n📦 正在将分片文件打包至 {FINAL_ZIP}...")
    with zipfile.ZipFile(FINAL_ZIP, 'w', zipfile.ZIP_DEFLATED) as z:
        for root, dirs, files in os.walk(OUTPUT_DIR):
            for file in files:
                z.write(os.path.join(root, file), file)
                # 可选：打包后删除原始分片以节省空间
                # os.remove(os.path.join(root, file))

    print(f"✨ 处理完成！总计下载记录数: {total_downloaded}")
    print(f"结果已保存至: {FINAL_ZIP}")

if __name__ == "__main__":
    main() # 可选：传入特定靶点名称，如 main("EGFR")

