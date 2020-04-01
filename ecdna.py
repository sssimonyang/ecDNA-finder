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

import pysam as ps

import utils


class EcDNA:
    def __init__(self, bam_file, split_read_mates, discordant_mates):
        self.bam = ps.AlignmentFile(bam_file, 'rb')
        self.split_read_mates = split_read_mates
        self.discordant_mates = discordant_mates
        self.max_insert = utils.extension
        self.read_length = 101
        self.long_intervals = []
        self.coverage_min = 10
        self.extend_size = 5000
        self.average_coverage = 5.20273
        self.coverage_std = 17.6794

        if self.split_read_mates and self.discordant_mates:
            self.build_long_intervals()

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
                if intersect_interval and self.similar_coverage(current_interval, intersect_interval):
                    whole_interval = self.combine(current_interval, intersect_interval)
                    if whole_interval:
                        current_support_mates.append(intersect_interval.mate)
                        long_intervals.remove(intersect_interval)
                        if intersect_interval == other_long_interval:
                            extend = False
                        else:
                            current_interval = intersect_interval.other_long_interval
                            long_intervals.remove(intersect_interval.other_long_interval)
                        current_result.append(whole_interval)
                        continue
                if self.extend_amplified(extend_interval):
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
            for result, support_mate in results:
                if len(result) == len(support_mate):
                    circ_results.append((result, support_mate))

    def similar_coverage(self, interval1, interval2):
        return True

    def extend_amplified(self, interval):
        return

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
            return utils.Interval(current_interval.chrom, current_interval.start, intersect_interval.end, strand='+')
        if current_interval.clipped == "right" and intersect_interval.clipped == 'left':
            return utils.Interval(current_interval.chrom, intersect_interval.start, current_interval.end, strand='-')
        return False

    def interval_coverage(self, interval):
        coverage = sum([sum(a) for a in self.bam.count_coverage(interval.chrom, interval.start, interval.end)]) / (
                interval.end - interval.start)
        return coverage

    def window_coverage(self, i, window_size=-1):
        if window_size == -1:
            window_size = self.max_insert - self.read_length

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
        self.long_intervals.sort(key=lambda i: i.raw_coverage, reverse=True)
