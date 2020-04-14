#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
@author: SSSimon Yang
@contact: yangjingkang@126.com
@file: main.py
@time: 2020/3/27 21:13
@desc:
"""

# import argparse
import multiprocessing as mp
import os
import pickle
import time

import extract
import mates
import utils
from ecdna import ECDNA

# parser = argparse.ArgumentParser(description="ecDNA-finder")
# parser.add_argument('-d', '--rundir', help='your run directory')
# args = parser.parse_args()

# prepare
# os.chdir(os.path.join("/public/home/zhangjing1/yangjk/ecDNA/data/", args.rundir))
if not os.path.exists('temp'):
    os.mkdir('temp')
sorted_query_name_bam_file = 'sorted_query_name.bam'
sorted_coordinate_bam_file = 'sorted_coordinate.bam'
extracted_file = 'extracted.bam'
utils.cores = mp.cpu_count()
start = time.localtime()
print(time.strftime("%Y %b %d %H:%M:%S", start))
print(f"use cores {utils.cores}")

# extract useful reads and calculate extension
utils.insert_mean, utils.insert_std, beds = extract.extract(sorted_query_name_bam_file, extracted_file,
                                                            sample_size=utils.sample_size,
                                                            bed_size=utils.bed_size, tlen_min=utils.tlen_min,
                                                            mapq_cutoff=utils.mapq_cutoff)
utils.extension = utils.insert_mean + utils.std * utils.insert_std

# reads extracted cluster and split
extracted_coordinate_file, peaks_file = extract.cluster(extracted_file, cluster_distance=utils.cluster_distance)
split_peaks = extract.split_interval(peaks_file, cutoff=utils.peak_value_cutoff, cores=utils.cores)
with open("temp/split_peaks.pickle", 'wb') as f:
    pickle.dump([utils.insert_mean, utils.insert_std, beds, split_peaks], f)
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
print(f"before process, {len(split_read_mates)} sr mates, {len(discordant_mates)} discordant mates")
with open("temp/before_mates.pickle", 'wb') as f:
    pickle.dump([split_read_mates, discordant_mates], f)
split_read_mates = [split_read_mate for split_read_mate in split_read_mates if
                    len(split_read_mate.support_split_reads) >= utils.split_read_cutoff]
split_read_mates = utils.process_all_sr_mates(split_read_mates)
split_read_mates, discordant_mates = utils.assign_discordant_mates(split_read_mates, discordant_mates)
split_read_mates = [split_read_mate for split_read_mate in split_read_mates if
                    len(split_read_mate.support_discordant_reads) >= utils.discordant_read_cutoff]
split_read_mates.sort(key=lambda i: len(i.support_split_reads) + len(i.support_discordant_reads), reverse=True)
discordant_mates.sort(key=lambda i: len(i.support_discordant_reads), reverse=True)
print(f"after process, {len(split_read_mates)} sr mates, {len(discordant_mates)} discordant mates")
# assemble
ecdna = ECDNA(sorted_coordinate_bam_file, split_read_mates, discordant_mates, beds, max_insert=utils.extension,
              extend_size=utils.extend_size)
with open("temp/ecdna_finder_temp.pickle", 'wb') as f:
    pickle.dump([ecdna.depth_median, ecdna.depth_average,
                 ecdna.depth_std, split_read_mates, discordant_mates], f)
circ_results, not_circ_results = ecdna.assemble()
utils.circ_result_out(circ_results)

# analyse

# write output
with open("temp/ecdna_finder_result.pickle", 'wb') as f:
    pickle.dump([circ_results, not_circ_results], f)
print('write success')
end = time.localtime()
print(time.strftime("%Y %b %d %H:%M:%S", end))
print(f'use time {((time.mktime(end) - time.mktime(start)) / 3600):.2f}h')
