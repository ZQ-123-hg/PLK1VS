#!/usr/bin/env python

import math
import argparse
import json
from glob import glob
from multiprocessing import Pool

import h5py
import numpy as np
import pandas as pd
from sklearn.externals import joblib

# 使用相对路径或环境变量
MODEL_BASE_PATH = './admetlab_models'
fp_list = ['MACCS', 'ECFP2_2048', 'ECFP4_1024', 'ECFP4_2048', 'ECFP6_2048']

usage = ('admetlab_pred.py --input_dir <dir_path> --output_dir <dir_path> [--threads <int> --suffix <str>]\n\n'
'Example: admetlab_pred.py -i ./input/ -o ./output/ -t 20 -sf 00')

parser = argparse.ArgumentParser(usage=usage)
parser.add_argument('-i', '--input_dir', required=True)
parser.add_argument('-o', '--output_dir', required=True)
parser.add_argument('-sf', '--suffix', type=str)
parser.add_argument('-t', '--threads', default=1, type=int)
args = parser.parse_args()


def lipinski(d):
    res = np.array([d['MW'] <= 500, d['LogP'] <= 5, d['naccr'] <= 10, d['ndonr'] <= 5])
    return res.mean(0)

def ghose(d):
    res = np.array([np.logical_and(d['LogP'] > 0.4, d['LogP'] < 5.6),
                    np.logical_and(d['MW'] > 160, d['MW'] < 480),
                    np.logical_and(d['MR'] > 40, d['MR'] < 130),
                    np.logical_and(d['nta'] > 20, d['nta'] < 70)])
    return res.mean(0)

def oprea(d):
    res = np.array([d['nring'] >= 3, d['nrigidbond'] >= 18, d['nRotbond'] >= 6])
    return res.mean(0)

def veber(d):
    res = np.array([d['nRotbond'] <= 10,
                    np.logical_or(d['TPSA'] <= 140, d['naccr'] + d['ndonr'] <= 12)])
    return res.mean(0)

def varma(d, r):
    res = np.array([d['MW'] <= 500,
                    d['TPSA'] <= 125,
                    np.logical_and(r['LogD7.4'] > 2, r['LogD7.4'] < 5),
                    d['naccr'] + d['ndonr'] <= 9,
                    d['nRotbond'] <= 12])
    return res.mean(0)


# 使用相对路径
df_meta = pd.read_csv(f'{MODEL_BASE_PATH}/admetlab_models.txt', delimiter='\t', index_col=0)
df_meta = df_meta[df_meta['Model path'].notnull()]
des_path = '/des_' + args.suffix + '.csv' if args.suffix else '/des.csv'
df_des = pd.read_csv(args.input_dir + des_path, index_col=0)
with open(f'{MODEL_BASE_PATH}/admetlab_labels.json', 'r') as f:
    label_dict = json.load(f)

fp_dict = {}
name_dict = {}
for fp_name in fp_list:
    fp_file_name = '/' + fp_name + '_' + args.suffix + '.h5' if args.suffix else '/' + fp_name + '.h5'
    f = h5py.File(args.input_dir + fp_file_name, 'r')
    fp_dict[fp_name] = np.array(f['array'])
    name_dict[fp_name] = list(np.array(f['name']))
    f.close()

df_des_name = list(df_des.index)
for fp_name, name_list in name_dict.items():
    # 验证名称一致性
    if name_list != df_des_name:
        print(f"Warning: Name mismatch for {fp_name}")
        print(f"Fingerprint names: {len(name_list)}, Descriptor names: {len(df_des_name)}")
        
        # 取交集
        common_names = set(name_list) & set(df_des_name)
        if len(common_names) == 0:
            raise ValueError("No common names found between fingerprints and descriptors")
        
        # 重新对齐数据
        name_order = [name for name in df_des_name if name in common_names]
        df_des = df_des.loc[name_order]
        for fn in fp_list:
            if fn in fp_dict:
                name_to_index = {name: i for i, name in enumerate(name_dict[fn])}
                indices = [name_to_index[name] for name in name_order if name in name_to_index]
                fp_dict[fn] = fp_dict[fn][indices]
        
        df_des_name = name_order
        break

split_len = int(math.ceil(float(len(df_des_name)) / args.threads))
index_list = [(i * split_len, min((i + 1) * split_len, len(df_des_name))) for i in range(args.threads)]


def predict(index):
    res = {}
    start, end = index
    prop_list = df_meta.index
    
    for prop in prop_list:
        task_type = df_meta.loc[prop, 'Type']
        model_path = df_meta.loc[prop, 'Model path']
        model_type = df_meta.loc[prop, 'Method']
        des_name = df_meta.loc[prop, 'Descriptor']
        print(f"Predicting {prop} using {model_type} with {des_name}")

        if des_name == '2D':
            des_array = np.array(df_des[label_dict[prop]][start:end])
            inf_index = np.where(np.sum(np.isfinite(des_array), 1) != des_array.shape[1])[0]
            des_array[inf_index] = 0
        else:
            des_array = fp_dict[des_name][start:end]

        if task_type == 'Regression':
            cf = joblib.load(f'{MODEL_BASE_PATH}/{model_path}')
            cf.verbose = 0  # 减少输出
            res[prop] = cf.predict(des_array)
            res[prop][inf_index] = np.nan
        elif task_type == 'Classification':
            model_paths = glob(f'{MODEL_BASE_PATH}/{model_path}')
            cf_list = [joblib.load(p) for p in model_paths]
            for cf in cf_list:
                cf.verbose = 0  # 减少输出
            res_list = [cf.predict_proba(des_array)[:, 1] for cf in cf_list]
            res[prop] = np.array(res_list).mean(0)
        else:
            raise Exception('Incorrect task type: ' + task_type)
        
    # 后处理计算
    res['LogP'] = df_des['LogP'][start:end]
    res['LD50 of acute toxicity'] = 10 ** -res['LD50 of acute toxicity'] * df_des['MW'][start:end] * 1000
    res['Lipinski'] = lipinski(df_des[start:end])
    res['Ghose'] = ghose(df_des[start:end])
    res['Oprea'] = oprea(df_des[start:end])
    res['Veber'] = veber(df_des[start:end])
    res['Varma'] = varma(df_des[start:end], res)
    
    return res


if __name__ == '__main__':
    pool = Pool(args.threads)
    res_list = pool.map(predict, index_list)
    pool.close()

    df_list = [pd.DataFrame(res) for res in res_list]
    df = pd.concat(df_list)
    df.index = df_des.index
    res_path = '/results_' + args.suffix + '.csv' if args.suffix else '/results.csv'
    df.to_csv(args.output_dir + res_path)
    print(f"Prediction completed. Results saved to {args.output_dir + res_path}")