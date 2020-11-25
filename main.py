#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# MIT License

# Copyright (c) 2020 SSSimon Yang

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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
parser.add_argument('-coord', '--coordinate', help='bam file sorted by coordinate', required=True)
parser.add_argument('-query', '--queryname', help='bam file sorted by queryname', required=True)
parser.add_argument('-dir', '--dirname', help='result directory', required=True)
parser.add_argument('-cutoff', '--cutoff', help='seed interval cutoff and support read cutoff', required=False,
                    default=0, type=int)
args = parser.parse_args()

# prepare
print(f'Start {args.dirname}')
if not os.path.exists(args.dirname):
    os.mkdir(args.dirname)
os.chdir(args.dirname)
sorted_query_name_bam_file = args.queryname
sorted_coordinate_bam_file = args.coordinate
utils.cutoff = args.cutoff
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
if utils.cutoff == 0:
    utils.cutoff = round(utils.depth_average / 20) or 1
print(f"cutoff {utils.cutoff}")
utils.extension = utils.insert_mean + utils.std * utils.insert_std
extracted_coordinate_file, peaks_file = extract.cluster(extracted_file, cluster_distance=utils.cluster_distance)
split_peaks = extract.split_interval(peaks_file, cutoff=utils.cutoff, cores=utils.cores)

with open("First.pickle", 'wb') as f:
    pickle.dump([extracted_coordinate_file, utils.extension, utils.depth_average, utils.depth_std, split_peaks], f)
end = time.localtime()
print(f'Now use time {((time.mktime(end) - time.mktime(start)) / 3600):.2f}h')

# find split read mates and discordant read mates
manager = mp.Manager()
pool = mp.Pool(processes=utils.cores)
split_read_mates = manager.list()
discordant_mates = manager.list()
mate = mates.Mates(extracted_coordinate_file, bam_file=sorted_coordinate_bam_file, extension=utils.extension,
                   mapq_cutoff=utils.mapq_cutoff,
                   interval_p_cutoff=utils.interval_p_cutoff)
for peaks in split_peaks:
    pool.apply_async(mate.find_mates, args=(peaks, split_read_mates, discordant_mates))
pool.close()
pool.join()
split_read_mates = list(split_read_mates)
discordant_mates = list(discordant_mates)

with open("Second.pickle", 'wb') as f:
    pickle.dump([sorted_coordinate_bam_file, split_read_mates, discordant_mates, utils.depth_average, utils.depth_std,
                 utils.extension], f)
end = time.localtime()
print(f'Now use time {((time.mktime(end) - time.mktime(start)) / 3600):.2f}h')

# filter and assemble
ecdna = ECDNA(sorted_coordinate_bam_file, split_read_mates, depth_average=utils.depth_average,
              depth_std=utils.depth_std, max_insert=utils.extension, extend_size=utils.extend_size,
              cutoff=utils.cutoff)
circ_results, not_circ_results = ecdna.assemble()
utils.circ_result_out(circ_results)
utils.not_circ_result_out(not_circ_results)

# analyse

# write output
with open("results.pickle", 'wb') as f:
    pickle.dump([circ_results, not_circ_results], f)
end = time.localtime()
print(time.strftime("%Y %b %d %H:%M:%S", end))
print(f'All use time {((time.mktime(end) - time.mktime(start)) / 3600):.2f}h')
