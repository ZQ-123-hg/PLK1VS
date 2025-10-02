import os
import pyarrow.parquet as pq
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.QED import qed
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
import sascorer
from concurrent.futures import ProcessPoolExecutor


# 读取文件并解析 SMILES 和化合物编号，同时跳过无效的 SMILES
def read_smiles_file(file_path):
    df = pd.read_csv(file_path, sep=' ', header=None, names=['SMILES', 'Compound_ID'])
    valid_smiles_list = []
    valid_compound_ids = []
    for smiles, compound_id in zip(df['SMILES'], df['Compound_ID']):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                valid_smiles_list.append(smiles)
                valid_compound_ids.append(compound_id)
            else:
                print(f"无法解析 SMILES 字符串: {smiles}, Compound_ID: {compound_id}")
        except Exception as e:
            print(f"处理 SMILES 字符串 {smiles} (Compound_ID: {compound_id}) 时出现错误: {e}")
    return valid_smiles_list, valid_compound_ids


# 创建结构警示过滤器
params = FilterCatalogParams()
params.AddCatalog(FilterCatalogParams.FilterCatalogs.ALL)
catalog = FilterCatalog(params)


# 定义一个函数来计算单个化合物的属性
def calculate_properties(smiles, compound_id):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            # 计算各属性
            MW = Descriptors.MolWt(mol)
            LogP = Descriptors.MolLogP(mol)
            HBA = rdMolDescriptors.CalcNumHBA(mol)
            HBD = rdMolDescriptors.CalcNumHBD(mol)
            RotatableBonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
            TPSA = rdMolDescriptors.CalcTPSA(mol)
            QED = qed(mol)
            AROMs = rdMolDescriptors.CalcNumAromaticRings(mol)
            ALERTS = catalog.HasMatch(mol)
            SAscore = sascorer.calculateScore(mol)
            StereoCenters = rdMolDescriptors.CalcNumAtomStereoCenters(mol)

            # 返回结果字典
            result = {
                'Compound_ID': compound_id,
                'SMILES': smiles,
                'MW': MW,
                'LogP': LogP,
                'HBA': HBA,
                'HBD': HBD,
                'RotatableBonds': RotatableBonds,
                'TPSA': TPSA,
                'QED': QED,
                'AROMs': AROMs,
                'ALERTS': ALERTS,
                'SAscore': SAscore,
                'StereoCenters': StereoCenters
            }
            return result
        else:
            print(f"无法解析 SMILES 字符串: {smiles}, Compound_ID: {compound_id}")
            return None
    except Exception as e:
        print(f"处理 SMILES 字符串 {smiles} (Compound_ID: {compound_id}) 时出现错误: {e}")
        return None


# 处理数据块
def process_chunk(chunk, chunk_index, file_prefix):
    # 直接解包 chunk
    smiles_list, compound_ids = chunk
    # 根据节点 CPU 核心数设置进程池大小
    ncpu = 96
    # 使用 85% 的核心
    max_workers = int(ncpu * 0.85)
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(calculate_properties, smiles, compound_id)
                   for smiles, compound_id in zip(smiles_list, compound_ids)]
        results = [future.result() for future in futures if future.result() is not None]
    # 将结果保存为 Parquet 文件
    result_df = pd.DataFrame(results)
    result_df.to_parquet(f'{file_prefix}_chunk_{chunk_index}.parquet', index=False)
    return f'{file_prefix}_chunk_{chunk_index}.parquet'


# 处理单个文件的主函数
def process_single_file(file_path):
    smiles_list, compound_ids = read_smiles_file(file_path)

    # 根据节点内存情况调整分块大小
    totalRAM = 251 * 1024  # 转换为 MB
    # 假设每个化合物处理大约需要 0.1 MB 内存
    memory_per_compound = 0.1
    # 每个进程预计占用的内存
    memory_per_process = memory_per_compound * 1000
    # 结合 CPU 核心使用数量和总内存确定分块大小
    ncpu = 96
    max_workers = int(ncpu * 0.85)
    chunk_size = int((totalRAM / memory_per_process) * max_workers)

    file_name = os.path.splitext(os.path.basename(file_path))[0]
    chunks = [(smiles_list[i:i + chunk_size], compound_ids[i:i + chunk_size])
              for i in range(0, len(smiles_list), chunk_size)]

    parquet_files = []
    for i, chunk in enumerate(chunks):
        parquet_file = process_chunk(chunk, i, file_name)
        parquet_files.append(parquet_file)

    # 检查文件是否存在
    for file in parquet_files:
        if not os.path.exists(file):
            print(f"文件 {file} 不存在，请检查。")
            break

    # 定义 CSV 文件路径
    csv_file = f'{file_name}.csv'

    # 分块处理每个 Parquet 文件
    chunk_size_parquet = 10000  # 可以根据实际情况调整分块大小

    # 处理第一个文件，写入表头
    first_table = pq.read_table(parquet_files[0])
    first_batches = first_table.to_batches(chunk_size_parquet)
    for i, batch in enumerate(first_batches):
        first_df = batch.to_pandas()
        if i == 0:
            first_df.to_csv(csv_file, index=False)
        else:
            first_df.to_csv(csv_file, mode='a', header=False, index=False)

    # 处理剩余的文件，追加数据
    for file in parquet_files[1:]:
        table = pq.read_table(file)
        batches = table.to_batches(chunk_size_parquet)
        for batch in batches:
            df = batch.to_pandas()
            df.to_csv(csv_file, mode='a', header=False, index=False)

    print(f"合并并转换完成，结果保存为 {csv_file}")

    # 删除中间的 Parquet 文件
    for file in parquet_files:
        os.remove(file)


if __name__ == '__main__':
    # 定义多个文件路径
    file_paths = ['test.txt']  # 替换为实际的文件路径

    for file_path in file_paths:
        print(f"开始处理文件: {file_path}")
        process_single_file(file_path)
        print(f"文件 {file_path} 处理完成。")

