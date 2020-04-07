#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
@author: SSSimon Yang
@contact: yangjingkang@126.com
@file: ecdna.py
@time: 2020/3/25 19:58
@desc:
"""
import copy

import numpy as np
import pysam as ps

import utils


class ECDNA:
    def __init__(self, bam_file, split_read_mates, discordant_mates, beds, max_insert=utils.extension, extend_size=5000,
                 std=2):
        self.bam = ps.AlignmentFile(bam_file, 'rb')
        self.split_read_mates = split_read_mates
        self.discordant_mates = discordant_mates
        self.beds = beds
        self.max_insert = max_insert
        self.extend_size = extend_size
        self.std = std
        self.long_intervals = []
        self.build_long_intervals()

        self.coverages = []
        for i in self.beds:
            self.coverages.append(self.interval_coverage(i))
        median_coverage = np.median(self.coverages)
        self.filter_coverage = [c for c in self.coverages if c < 2 * median_coverage]
        self.coverage_median = np.median(self.filter_coverage)  # 5.48
        self.coverage_average = np.average(self.filter_coverage)  # 5.59
        self.coverage_std = np.std(self.filter_coverage)  # 1.75
        self.upper_limit = self.coverage_average + self.coverage_std
        self.filter()

    def assemble(self):
        results = []
        long_intervals = copy.copy(self.long_intervals)
        while long_intervals:
            long_interval = long_intervals.pop(0)
            current_result = []
            current_support_mates = []
            other_long_interval = long_interval.other_long_interval
            current_interval = long_interval
            extend = True
            while extend:
                extend_interval = self.interval_extend(current_interval, self.extend_size)
                intersect_interval = extend_interval.intersects(long_intervals)
                if intersect_interval and intersect_interval.clipped == utils.opposite[current_interval.clipped]:
                    if intersect_interval.clipped == 'left':
                        extend_interval.start = intersect_interval.start
                    else:
                        extend_interval.end = intersect_interval.end
                    current_interval.extend_coverage.append(self.interval_coverage(extend_interval))
                    if self.similar_coverage(current_interval, intersect_interval):
                        whole_interval = self.combine(current_interval, intersect_interval)
                        current_support_mates.append(intersect_interval.mate)
                        long_intervals.remove(intersect_interval)
                        if intersect_interval == other_long_interval:
                            extend = False
                        else:
                            current_interval = intersect_interval.other_long_interval
                            long_intervals.remove(intersect_interval.other_long_interval)
                        current_result.append(whole_interval)
                        continue
                if self.extend_amplified(current_interval, extend_interval):
                    if current_interval.clipped == 'left':
                        current_interval.end = extend_interval.end
                    else:
                        current_interval.start = extend_interval.start
                    current_interval.extend_coverage.append(self.interval_coverage(extend_interval))
                    current_interval.coverage = self.interval_coverage(current_interval)
                else:
                    current_result.append(current_interval)
                    extend = False
            results.append((current_result, current_support_mates))

        circ_results = []
        not_circ_results = []
        for result, support_mate in results:
            if len(result) == len(support_mate):
                circ_results.append((result, support_mate))
            else:
                not_circ_results.append((result, support_mate))
        return circ_results, not_circ_results

    def similar_coverage(self, interval1, interval2):
        extend_coverage = np.array(interval1.extend_coverage)
        low = np.median(extend_coverage) - self.std * self.coverage_std
        high = np.median(extend_coverage) + self.std * self.coverage_std
        if np.all(low <= extend_coverage) and np.all(extend_coverage <= high) and low <= interval2.raw_coverage <= high:
            return True
        else:
            return False

    def extend_amplified(self, current_interval, extend_interval):
        extend_coverage = self.interval_coverage(extend_interval)
        if extend_coverage < self.upper_limit:
            return False
        median = np.median(current_interval.extend_coverage)
        low = median - self.std * self.coverage_std
        high = median + self.std * self.coverage_std
        if low <= extend_coverage <= high and low <= current_interval.raw_coverage <= high:
            return True
        else:
            return False

    @staticmethod
    def interval_extend(long_interval, extend_size):
        if long_interval.clipped == 'left':
            return utils.Interval(long_interval.chrom, long_interval.end,
                                  long_interval.end + extend_size)
        if long_interval.clipped == 'right':
            return utils.Interval(long_interval.chrom, long_interval.start - extend_size,
                                  long_interval.start)

    @staticmethod
    def combine(current_interval, intersect_interval):
        if current_interval.clipped == "left" and intersect_interval.clipped == 'right':
            interval = utils.WholeInterval(current_interval.chrom, current_interval.start, intersect_interval.end,
                                           strand='+', coverages=current_interval.extend_coverage)
        if current_interval.clipped == "right" and intersect_interval.clipped == 'left':
            interval = utils.WholeInterval(current_interval.chrom, intersect_interval.start, current_interval.end,
                                           strand='-', coverages=current_interval.extend_coverage)
        return interval

    def interval_coverage(self, interval):
        coverage = sum([sum(a) for a in self.bam.count_coverage(interval.chrom, interval.start, interval.end)]) / (
                interval.end - interval.start)
        return coverage

    def window_coverage(self, i, window_size=-1):
        if window_size == -1:
            window_size = self.max_insert

        def win_breakup(i, window_size):
            for k in range(i.start, i.end, window_size):
                yield utils.Interval(i.chrom, k, k + window_size - 1)

        for k in win_breakup(i, window_size):
            yield (k, self.interval_coverage(k))

    def build_long_intervals(self):
        long_interval1s = [utils.LongInterval(i.interval1, i.clipped1, i) for i in self.split_read_mates]
        long_interval2s = [utils.LongInterval(i.interval2, i.clipped2, i) for i in self.split_read_mates]
        for long_interval1, long_interval2 in zip(long_interval1s, long_interval2s):
            long_interval1.other_long_interval = long_interval2
            long_interval2.other_long_interval = long_interval1
        self.long_intervals.extend(long_interval1s)
        self.long_intervals.extend(long_interval2s)
        for long_interval in self.long_intervals:
            long_interval.raw_coverage = self.interval_coverage(long_interval)
            long_interval.extend_coverage.append(long_interval.raw_coverage)
        self.long_intervals.sort(key=lambda i: i.raw_coverage, reverse=True)

    def filter(self):
        self.long_intervals = [long_interval for long_interval in self.long_intervals if
                               long_interval.raw_coverage > self.upper_limit]


# todo calculate score

if __name__ == '__main__':
    import pickle

    with open("/public/home/zhangjing1/yangjk/ecDNA/data/ecdna_3427/results.pickle", 'rb') as f:
        split_read_mates, discordant_mates = pickle.load(f)
    bam_file = '/public/home/zhangjing1/yangjk/ecDNA/data/sorted_coordinate.bam'
    ecDNA = EcDNA(bam_file, split_read_mates, discordant_mates)
    results = ecDNA.assemble()
    print(results)
