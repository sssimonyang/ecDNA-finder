#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
@author: SSSimon Yang
@contact: yangjingkang@126.com
@file: analysis.py
@time: 2020/5/14 12:48
@desc:
"""
import os

import pandas as pd

path = "/public/home/zhangjing1/yangjk/ecDNA/result"
os.chdir(path)
dirs = os.listdir(path)
circ_number = 0
result = pd.DataFrame()
for i in dirs:
    if not i.startswith('F'):
        continue
    if not os.path.exists(os.path.join(path, i, 'circ_results.tsv')):
        continue
    data = pd.read_csv(os.path.join(path, i, 'circ_results.tsv'), sep='\t')
    if data.empty:
        continue
    else:
        print(i)
        data.circ_id += circ_number
        data['source'] = i
        result = pd.concat([result, data])
    circ_number = result.circ_id.max()
result.reset_index(drop=True, inplace=True)
result.to_csv('f_all_circ_result.tsv', sep='\t', index=False)
