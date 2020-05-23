#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: SSSimon Yang
@contact: yangjingkang@126.com
@file: mates.py
@time: 2020/2/16 12:37
@desc:
"""

import pandas as pd
import pysam as ps

import utils


class Mates:

    def __init__(self, extracted_bam_file, bam_file, extension, mapq_cutoff=40, interval_p_cutoff=0.5):
        self.extracted_bam_file = extracted_bam_file
        self.bam_file = bam_file
        self.extension = extension
        self.mapq_cutoff = mapq_cutoff
        self.interval_p_cutoff = interval_p_cutoff

    def interval_depth(self, bam, interval):
        depth = sum([sum(a) for a in bam.count_coverage(interval.chrom, interval.start, interval.end)]) / (
                interval.end - interval.start)
        return depth

    def find_mates(self, peaks, split_read_mates, discordant_mates):
        peaks = pd.DataFrame.from_records(peaks, columns=['chrom', 'start', 'end'])
        extracted_bam = ps.AlignmentFile(self.extracted_bam_file, "rb")

        for index, interval in peaks.iterrows():
            interval = utils.to_interval(interval)
            candidate_mates = utils.get_mate_intervals(extracted_bam, interval, self.mapq_cutoff)
            if len(candidate_mates) == 0:
                continue
            realignment_interval_extended = utils.get_realignment_intervals(candidate_mates, self.interval_p_cutoff)
            if realignment_interval_extended.empty:
                continue
            for _, mate_interval in realignment_interval_extended.iterrows():
                mate_interval = utils.to_interval(mate_interval)
                iteration_split_read_mates = []
                support_discordant_reads = [(read, extracted_bam.mate(read)) for read in
                                            extracted_bam.fetch(interval.chrom, interval.start, interval.end)
                                            if mate_interval.contain(utils.to_interval(extracted_bam.mate(read)))]
                if len(support_discordant_reads) <= 1:
                    continue
                support_clipped_reads = [read for read in
                                         extracted_bam.fetch(interval.chrom, interval.start, interval.end)
                                         if utils.is_clipped(read)]
                support_sr_reads = [read for read in
                                    extracted_bam.fetch(interval.chrom, interval.start, interval.end)
                                    if utils.is_clipped(read) and read.has_tag('SA')]
                if support_sr_reads:
                    discordant = False
                else:
                    if support_clipped_reads:
                        remap_info = utils.get_longest_soft_clipped_bases(support_clipped_reads)
                        if remap_info[1] < 10:
                            discordant = True
                        else:
                            discordant = False
                    else:
                        discordant = True

                if discordant:
                    discordant_mate = utils.Discordant_Mate(interval, mate_interval)
                    discordant_mate.support_discordant_reads = utils.show_mates(support_discordant_reads, sr=False)
                    discordant_mates.append(discordant_mate)
                    continue

                if support_sr_reads:
                    for read in support_sr_reads:
                        supplementary = utils.SA(read.get_tag("SA"))
                        if supplementary.mapping_quality > self.mapq_cutoff and mate_interval.contain(
                                supplementary):
                            if utils.clipped_direction(supplementary) == 'right':
                                clipped_supplementary = 'right'
                                interval_supplementary = utils.Interval(supplementary.chrom,
                                                                        min(supplementary.end - self.extension,
                                                                            mate_interval.start),
                                                                        supplementary.end,
                                                                        supplementary.strand)
                            else:
                                clipped_supplementary = 'left'
                                interval_supplementary = utils.Interval(supplementary.chrom,
                                                                        supplementary.start,
                                                                        max(supplementary.start + self.extension,
                                                                            mate_interval.end),
                                                                        supplementary.strand)
                            if utils.clipped_direction(read) == 'right':
                                clipped_primary = 'right'
                                interval_primary = utils.Interval(read.reference_name,
                                                                  read.reference_end - self.extension,
                                                                  read.reference_end,
                                                                  '-' if read.is_reverse else '+')
                            else:
                                clipped_primary = 'left'
                                interval_primary = utils.Interval(read.reference_name,
                                                                  read.reference_start,
                                                                  read.reference_start + self.extension,
                                                                  '-' if read.is_reverse else '+')
                            split_read_mate = utils.SR_Mate(interval_primary, clipped_primary,
                                                            interval_supplementary, clipped_supplementary)
                            split_read_mate.from_read = utils.show_mates([(read, supplementary)], sr=True)[0]
                            iteration_split_read_mates.append(split_read_mate)
                else:
                    if support_clipped_reads:
                        remap_info = utils.get_longest_soft_clipped_bases(support_clipped_reads)
                        if remap_info[1] >= 10:
                            pass
                            # print(f'{remap_info[2].query_name} need remap {interval} {mate_interval}')
                            # todo edlib 利用remap_info remap获得sr_mate both = True时考虑双向
                if iteration_split_read_mates:
                    iteration_split_read_mates = utils.process_sr_mate(iteration_split_read_mates)

                # second pass to add discordant read info
                if iteration_split_read_mates:
                    iteration_split_read_mates = utils.assign_discordant_reads(iteration_split_read_mates,
                                                                               support_discordant_reads)
                    bam = ps.AlignmentFile(self.bam_file, 'rb')
                    for split_read_mate in iteration_split_read_mates:
                        split_read_mate.interval1.depth = self.interval_depth(bam, split_read_mate.interval1)
                        split_read_mate.interval2.depth = self.interval_depth(bam, split_read_mate.interval2)
                    split_read_mates.extend(iteration_split_read_mates)

        extracted_bam.close()
        return
