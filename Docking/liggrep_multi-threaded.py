import sys
import concurrent.futures
import glob
from multiprocessing import cpu_count

from liggrep import main as liggrep_main

def process_ligands(receptor, ligand, filters):
    print(f"Processing ligand: {ligand}")
    liggrep_main([receptor, ligand, filters])
    print(f"Processed ligand: {ligand}")

if __name__ == "__main__":
    # 检查命令行参数
    if len(sys.argv) != 4:
        print("Usage: python your_script.py receptor ligands filters")
        sys.exit(1)

    # 读取命令行参数
    receptor, ligands_pattern, filters = sys.argv[1:]

    # 获取匹配的 ligands 文件
    ligands_files = glob.glob(ligands_pattern)

    # 指定并发执行的最大进程数
    num_processes = min(cpu_count(), 60)  # 使用最大可用CPU核心数，最多不超过40个进程

    # 创建进程池
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_processes) as executor:
        # 提交子任务给进程池
        # 每个子任务处理一个 ligands 文件
        futures = [executor.submit(process_ligands, receptor, ligand, filters) for ligand in ligands_files]

        # 等待所有子任务完成
        concurrent.futures.wait(futures)
