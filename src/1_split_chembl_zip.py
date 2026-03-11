import os
import zipfile
import pandas as pd
import io

import multiprocessing as mp

def process_chembl_zip_to_batches(zip_path, output_dir, chunk_size):
        """
        专门针对 ChEMBL 等大型压缩包的流式处理函数
        """
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        with zipfile.ZipFile(zip_path, 'r') as z:
            # 1. 获取压缩包内文件名（通常只有一个主数据文件）
            all_files = z.namelist()
            target_file = [f for f in all_files if f.endswith('.csv') or f.endswith('.tsv')][0]

            print(f"检测到目标文件: {target_file}")

            # 2. 开启流式读取，不解压到硬盘
            with z.open(target_file) as f:
                reader = pd.read_csv(
                    f, 
                    sep=None,
                    engine='python', 
                    chunksize=chunk_size, 
                    encoding='utf-8', 
                    on_bad_lines='warn',
                    low_memory=True
                )

                for i, chunk in enumerate(reader):
                    output_file = os.path.join(output_dir, f"chembl_part_{i+1}.csv")

                    # 3. 规范化保存：使用 utf-8-sig 确保 Excel/记事本打开不乱码
                    chunk.to_csv(
                        output_file, 
                        index=False, 
                        encoding='utf-8-sig', 
                        quoting=1
                    )

                    print(f"成功导出分片 {i+1}，包含 {len(chunk)} 行数据")


if __name__ == "__main__":

    process_chembl_zip_to_batches(
        zip_path='./datasets/chembl_36.zip', 
        output_dir='./chembl_split_data_zip',
        chunk_size=20000
    )