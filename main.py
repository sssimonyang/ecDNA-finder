#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
@author: SSSimon Yang
@contact: yangjingkang@126.com
@file: main.py
@time: 2020/3/27 21:13
@desc:
"""

import multiprocessing as mp
import os
import pickle
import time

import extract
import mates
import utils

utils.cores = mp.cpu_count()
os.chdir("/public/home/zhangjing1/yangjk/ecDNA/data/")
start = time.localtime()
print(time.strftime("%Y %b %d %H:%M:%S", start))
print(f"pid {utils.pid} file in directory ecdna_{utils.pid}")
print(f"use cores {utils.cores}")
os.mkdir(f'ecdna_{utils.pid}')
os.chdir(f'ecdna_{utils.pid}')
sorted_query_name_bam_file = '/public/home/zhangjing1/yangjk/ecDNA/data/sorted_query_name.bam'
extracted_file = 'extracted.bam'
mean, std = extract.extract(sorted_query_name_bam_file, extracted_file, sample_size=utils.sample_size,
                            tlen_min=utils.tlen_min, mapq_cutoff=utils.mapq_cutoff)
utils.insert_mean = mean
utils.insert_std = std
utils.extension = mean + utils.std * std

extracted_coordinate_file, peaks_file = extract.cluster(extracted_file, cluster_distance=utils.cluster_distance)
split_peaks = extract.split_interval(peaks_file, cores=utils.cores)
mate = mates.Mates(extracted_coordinate_file, extension=utils.extension, mapq_cutoff=utils.mapq_cutoff,
                   interval_p_cutoff=utils.interval_p_cutoff)

manager = mp.Manager()
pool = mp.Pool(processes=utils.cores)
split_read_mates = manager.list()
discordant_mates = manager.list()
for peaks in split_peaks:
    pool.apply_async(mate.find_mates, args=(peaks, split_read_mates, discordant_mates))
pool.close()
pool.join()
split_read_mates = list(split_read_mates)
discordant_mates = list(discordant_mates)
split_read_mates, discordant_mates = utils.assign_discordant_mates(split_read_mates, discordant_mates)
split_read_mates.sort(key=lambda i: len(i.support_split_reads) + len(i.support_discordant_reads), reverse=True)
discordant_mates.sort(key=lambda i: len(i.support_discordant_reads), reverse=True)

print(f"after process, {len(split_read_mates)} sr mates, {len(discordant_mates)} discordant mates")
with open("results.pickle", 'wb') as f:
    pickle.dump([split_read_mates, discordant_mates], f)
print('write success')
end = time.localtime()
print(time.strftime("%Y %b %d %H:%M:%S", end))
print(f'use time {((time.mktime(end) - time.mktime(start)) / 3600):.2f}h')
