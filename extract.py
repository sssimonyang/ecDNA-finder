#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: SSSimon Yang
@contact: yangjingkang@126.com
@file: extract.py
@time: 2020/3/10 16:22
@desc:
"""

import subprocess as sp
from random import random

import numpy as np
import pandas as pd
import pysam as ps

import utils


# nohup samtools sort -n -@ 4 -o sorted_query_name.bam sorted_coordinate.bam &

def extract(in_bam_file, out_bam_file, sample_size, bed_size, tlen_min, mapq_cutoff):
    in_bam = ps.AlignmentFile(in_bam_file, 'rb')
    utils.chrom_names = in_bam.references[:24]
    utils.chrom_lengths = in_bam.lengths[:24]
    out_bam = ps.AlignmentFile(out_bam_file, 'wb', template=in_bam)

    def process_same_qname_reads(same_qname_reads):
        write = False
        nonlocal all_num, write_num, discordant_reads_num, split_reads_num, insert_length, sample_num, beds, bed_num
        all_num += 1
        if any([read.reference_name not in utils.chrom_names or read.is_unmapped or read.is_qcfail or read.is_duplicate
                for read in same_qname_reads]):
            return
        if any([read.is_read1 for read in same_qname_reads]) and any([read.is_read2 for read in same_qname_reads]):
            if len(same_qname_reads) == 2:
                if all([read.mapq >= mapq_cutoff for read in same_qname_reads]):
                    if same_qname_reads[0].reference_name != same_qname_reads[1].reference_name or abs(
                            same_qname_reads[0].tlen) > tlen_min:
                        discordant_reads_num += 1
                        write = True
                    else:
                        if sample_num < sample_size and random() > 0.1:
                            insert_length.append(abs(same_qname_reads[0].tlen))
                            sample_num += 1
                        if bed_num < bed_size and random() > 0.1:
                            beds.append(utils.to_interval(same_qname_reads[0]))
                            bed_num += 1

            elif len(same_qname_reads) == 3:
                read_1 = [read for read in same_qname_reads if read.is_read1]
                read_2 = [read for read in same_qname_reads if read.is_read2]
                if len(read_1) == 2:
                    read = [read for read in read_1 if not read.is_supplementary][0]
                    supplementary = [read for read in read_1 if read.is_supplementary][0]
                    read_mate = read_2[0]
                else:
                    read = [read for read in read_2 if not read.is_supplementary][0]
                    supplementary = [read for read in read_2 if read.is_supplementary][0]
                    read_mate = read_1[0]
                if read.mapping_quality > mapq_cutoff and read_mate.mapping_quality > mapq_cutoff:
                    if utils.colinear(read, read_mate, supplementary):
                        split_reads_num += 1
                        write = True
            else:
                pass
        if write:
            write_num += 1
            for read in same_qname_reads:
                out_bam.write(read)

    all_num = write_num = discordant_reads_num = split_reads_num = sample_num = bed_num = 0
    beds = []
    insert_length = []
    temp_read = next(in_bam)
    query_name = temp_read.query_name
    same_qname_reads = [temp_read]
    for temp_read in in_bam:
        if temp_read.query_name != query_name:
            process_same_qname_reads(same_qname_reads)
            query_name = temp_read.query_name
            same_qname_reads = [temp_read]
        else:
            same_qname_reads.append(temp_read)
    process_same_qname_reads(same_qname_reads)
    in_bam.close()
    out_bam.close()
    mean = np.mean(insert_length)
    std = np.std(insert_length)

    print(f"all {all_num}\nextract {write_num}\n"
          f"discordant_reads {discordant_reads_num}\nsplit_reads {split_reads_num}\n"
          f"insert length mean {mean} std {std}")
    return mean, std, beds


def cluster(extracted_file, cluster_distance):
    # samtools sort -n -o sorted_query_name_extracted.bam extracted.bam
    # samtools sort -o sorted_coordinate_name_extracted.bam extracted.bam
    sp.call(f"samtools sort -n -o sorted_query_name_{extracted_file} {extracted_file}", shell=True)
    sp.call(f"samtools sort -o sorted_coordinate_{extracted_file} {extracted_file}", shell=True)
    sp.call(f"samtools index sorted_coordinate_{extracted_file}", shell=True)
    sp.call(
        f"bedtools genomecov -bg -ibam sorted_coordinate_{extracted_file} | sort -k 1,1 -k2,2n | "
        f"mergeBed -d {cluster_distance} -c 4 -o mean | sort -r -n -k 4,4 > peaks.bed", shell=True)

    return f"sorted_coordinate_{extracted_file}", 'peaks.bed'


def split_interval(peaks_file, cores):
    peaks = pd.read_csv(peaks_file, sep='\t', names=['chrom', 'start', 'end', 'value'],
                        dtype={'chrom': 'str', 'start': 'int', 'end': 'int', 'value': 'float'})
    print(f"interval {peaks.shape[0]}")
    chunks = cores * 100
    counter = 0
    split_peaks = []
    for i in range(0, chunks):
        split_peaks.append([])

    for _, interval in peaks.iterrows():
        if counter == chunks:
            counter = 0
        interval.start = interval.start - 100
        interval.end = interval.end + 100
        if interval.end - interval.start > 500:
            w_start = interval.start
            while w_start < interval.end:
                splitted = [interval.chrom, w_start, w_start + 300]
                w_start += 300
                if counter == chunks:
                    counter = 0
                    split_peaks[counter].append(splitted)
                else:
                    split_peaks[counter].append(splitted)
                    counter += 1
        else:
            split_peaks[counter].append([interval.chrom, interval.start, interval.end])
            counter += 1
    print(f"after split {len(split_peaks)} chunks {sum([len(i) for i in split_peaks])} intervals")
    return split_peaks
