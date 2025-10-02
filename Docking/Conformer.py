# -*- coding: utf-8 -*-
import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.transform import Rotation as R
from openbabel import openbabel

# 设置 Open Babel 的路径
OBABEL_PATH = '/path/obabel'

# 读取SMILES文件并生成3D结构
def smiles_to_3d(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.UFFOptimizeMolecule(mol)
    return mol

# 计算分子的中心点（几何中心）
def get_center_of_mass(mol):
    conf = mol.GetConformer()
    positions = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    center = np.mean(positions, axis=0)
    return center

# 绕中心点旋转分子
def rotate_molecule(mol, center, rotation_matrix):
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        pos = np.array(conf.GetAtomPosition(i)) - center
        new_pos = np.dot(rotation_matrix, pos) + center
        conf.SetAtomPosition(i, new_pos)

# 将RDKit分子对象转换为Open Babel分子对象
def rdkit_to_openbabel(mol):
    block = Chem.MolToMolBlock(mol)
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat("mol")
    obMol = openbabel.OBMol()
    if not obConversion.ReadString(obMol, block):
        raise ValueError("Failed to convert RDKit molecule to Open Babel")
    return obMol

# 保存为PDBQT文件
def save_as_pdbqt(obmol, filename):
    obConversion = openbabel.OBConversion()
    obConversion.SetOutFormat("pdbqt")
    if not obConversion.WriteFile(obmol, filename):
        raise IOError(f"Failed to write PDBQT file: {filename}")

# 主函数：读取SMILES码文件，生成3D结构，旋转并保存
def process_smiles_file(smiles_file):
    base_output_dir = 'output/test'
    if not os.path.exists(base_output_dir):
        os.makedirs(base_output_dir)

    with open(smiles_file, 'r') as f:
        smiles_list = f.readlines()

    for idx, line in enumerate(smiles_list):
        line = line.strip()
        parts = line.split()  # 假设第二列是我们需要的前缀
        if len(parts) < 2:
            print(f"Skipping line {idx}: {line} (not enough columns)")
            continue
            
        smiles = parts[0]  # SMILES在第一列
        prefix = parts[1]  # 前缀在第二列

        try:
            mol = smiles_to_3d(smiles)
            center = get_center_of_mass(mol)

            # 将RDKit分子对象转换为Open Babel分子对象
            obmol = rdkit_to_openbabel(mol)

            # 保存原始的PDBQT文件
            original_pdbqt = os.path.join(base_output_dir, f"{prefix}_original.pdbqt")
            save_as_pdbqt(obmol, original_pdbqt)

            # 进行100次随机旋转并保存
            for i in range(100):
                rotation_matrix = R.random().as_matrix()
                rotate_molecule(mol, center, rotation_matrix)
                obmol = rdkit_to_openbabel(mol)  # 更新Open Babel分子对象
                rotated_pdbqt = os.path.join(base_output_dir, f"{prefix}_rotated_{i}.pdbqt")
                save_as_pdbqt(obmol, rotated_pdbqt)
        except Exception as e:
            print(f"Error processing SMILES {smiles}: {e}")

if __name__ == "__main__":
    # 设置Open Babel的环境变量
    os.environ['PATH'] += os.pathsep + '/path/bin'

    # 使用80个线程
    os.environ['OB_THREAD_COUNT'] = '80'

    # 假设你的SMILES文件名为'smiles.txt'
    process_smiles_file('/input/test.smi')