#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
@author: SSSimon Yang
@contact: yangjingkang@126.com
@file: main.py
@time: 2020/3/27 21:13
@desc:
"""

import argparse
import multiprocessing as mp
import os
import pickle
import time

import extract
import mates
import utils
from ecdna import ECDNA

parser = argparse.ArgumentParser(description="ecDNA-finder")
parser.add_argument('-c', '--coordinate', help='bam file sorted by coordinate', required=True)
parser.add_argument('-q', '--queryname', help='bam file sorted by queryname', required=True)
parser.add_argument('-d', '--dirname', help='result directory', required=True)
args = parser.parse_args()

# prepare
print(f'start {args.dirname}')
if not os.path.exists(args.dirname):
    os.mkdir(args.dirname)
os.chdir(args.dirname)
sorted_query_name_bam_file = args.queryname
sorted_coordinate_bam_file = args.coordinate
extracted_file = 'extracted.bam'
utils.cores = mp.cpu_count()
start = time.localtime()
print(time.strftime("%Y %b %d %H:%M:%S", start))
print(f"use cores {utils.cores}")

# extract, cluster and split
utils.insert_mean, utils.insert_std, beds = extract.extract(sorted_query_name_bam_file, extracted_file,
                                                            sample_size=utils.sample_size,
                                                            bed_size=utils.bed_size, tlen_min=utils.tlen_min,
                                                            mapq_cutoff=utils.mapq_cutoff)
utils.depth_average, utils.depth_std = extract.compute_depth(sorted_coordinate_bam_file, beds)
### 谨慎更改
utils.peak_value_cutoff = round(utils.depth_average / 10) or 1
utils.extension = utils.insert_mean + utils.std * utils.insert_std
extracted_coordinate_file, peaks_file = extract.cluster(extracted_file, cluster_distance=utils.cluster_distance)
split_peaks = extract.split_interval(peaks_file, cutoff=utils.peak_value_cutoff, cores=utils.cores)

# find split read mates and discordant read mates
manager = mp.Manager()
pool = mp.Pool(processes=utils.cores)
split_read_mates = manager.list()
discordant_mates = manager.list()
mate = mates.Mates(extracted_coordinate_file, extension=utils.extension, mapq_cutoff=utils.mapq_cutoff,
                   interval_p_cutoff=utils.interval_p_cutoff)
for peaks in split_peaks:
    pool.apply_async(mate.find_mates, args=(peaks, split_read_mates, discordant_mates))
pool.close()
pool.join()
split_read_mates = list(split_read_mates)
discordant_mates = list(discordant_mates)

# filter and assemble
ecdna = ECDNA(sorted_coordinate_bam_file, split_read_mates, discordant_mates, depth_average=utils.depth_average,
              depth_std=utils.depth_std, max_insert=utils.extension, extend_size=utils.extend_size)
with open("read_mates.pickle", 'wb') as f:
    pickle.dump([ecdna.split_read_mates, ecdna.discordant_mates], f)
circ_results, not_circ_results = ecdna.assemble()
with open("results.pickle", 'wb') as f:
    pickle.dump([circ_results, not_circ_results], f)
utils.circ_result_out(circ_results)
utils.not_circ_result_out(not_circ_results)

# analyse

# write output
with open("results.pickle", 'wb') as f:
    pickle.dump([circ_results, not_circ_results], f)
print('write success')
end = time.localtime()
print(time.strftime("%Y %b %d %H:%M:%S", end))
print(f'use time {((time.mktime(end) - time.mktime(start)) / 3600):.2f}h')
