# 导入依赖包
from rdkit.Chem import AllChem as ch
from rdkit.Chem import Draw as d
import glob
import os

# 载入substructure结构
substructure_smiles = [
    'CNC1=NC=CC=N1', 'C12=NC=NC=C1C=CC=C2', 'C12=CC=CC=C1C=CC=N2', 
    'C1(NC=C2)=C2C=NC=N1', 'CNC1=NC=CC=C1', 'O=C(N1)CC2=C1C=CC=C2', 
    'C12=NC=NC=C1N=CN2'
]

# 获取所有SMI文件路径（使用相对路径）
smi_files = glob.glob('./input/*.smi')  # 修改为相对路径

# 对每个SMI文件分别进行处理
for smi_file in smi_files:
    # 获取文件名（不含路径）
    file_basename = os.path.basename(smi_file).replace('.smi', '')
    
    # 读取当前SMI文件中的所有分子
    with open(smi_file, 'r') as f:
        lines = f.readlines()
        mols = [ch.MolFromSmiles(line.strip()) for line in lines if ch.MolFromSmiles(line.strip()) is not None]

    print(f"文件 {file_basename}.smi 读取了 {len(mols)} 个分子")

    # 对每个子结构分别进行搜索并保存匹配的结构
    for idx, smiles in enumerate(substructure_smiles):
        # 将SMILES转化为Mol对象
        pattern = ch.MolFromSmiles(smiles)
        if pattern is not None:
            print(f"子结构 {idx+1} 创建成功！")
            # 匹配包含当前子结构的分子
            matching_molecules = [mol for mol in mols if mol.HasSubstructMatch(pattern)]
            print(f"匹配子结构 {idx+1} 的分子数:", len(matching_molecules))

            # 保存匹配的结构到SD文件（使用相对路径）
            output_filename = f'./output/matched_structures_{file_basename}_substructure_{idx+1}.sdf'
            
            # 确保输出目录存在
            os.makedirs(os.path.dirname(output_filename), exist_ok=True)
            
            with ch.SDWriter(output_filename) as w:
                for mol in matching_molecules:
                    w.write(mol)
            print(f"结果已保存至: {output_filename}")
        else:
            print(f"子结构 {idx+1} 创建失败！")

print("所有处理完成！")