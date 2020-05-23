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
    def __init__(self, bam_file, split_read_mates, depth_average=utils.depth_average,
                 depth_std=utils.depth_std, max_insert=utils.extension, extend_size=utils.extend_size, cutoff=False):
        self.bam = ps.AlignmentFile(bam_file, 'rb')
        self.split_read_mates = split_read_mates
        self.depth_average = depth_average
        self.depth_std = depth_std
        self.max_insert = max_insert
        self.extend_size = extend_size
        self.upper_limit = round(self.depth_average + 1.5 * self.depth_std)
        self.cutoff = cutoff
        if self.cutoff:
            self.split_read_cutoff = round(self.depth_average / 20) or 1
            self.discordant_read_cutoff = round(self.depth_average / 20)
            print(f"split_read_cutoff {self.split_read_cutoff}, discordant_read_cutoff {self.discordant_read_cutoff}")
        self.filter()

        self.long_intervals_for_index = []
        self.mates = []
        self.build_long_intervals()

    def assemble(self):
        results = []
        while self.mates:
            long_interval1, long_interval2, first_or_not = self.mates.pop(0)
            circ_intervals = []
            support_mates = []
            currents = [(long_interval1, [long_interval1.other_index_interval], circ_intervals, support_mates, False)]
            circ_first = False
            while currents:
                current = currents.pop()
                for extend in self.extend(*current[:-1]):
                    if extend[0]:
                        currents.append(extend)
                    else:
                        results.append((extend[2], extend[3]))
                    if extend[-1] == 0:
                        circ_first = True
                # if currents:
                #     print(current[-1])
                # else:
                #     print([])
            if not circ_first and first_or_not:
                # todo
                self.mates.append((long_interval2, long_interval1, False))
            all_support_mates = [support_mate for result in results for support_mate in result[1]]
            self.mates = [mate for mate in self.mates if mate[0].mate not in all_support_mates]
            # print(f'left {len(self.mates)} sr mates')
        circ_results = []
        not_circ_results = []
        for circ_intervals, support_mates in results:
            if len(circ_intervals) == len(support_mates):
                circ_results.append((circ_intervals, support_mates))
            else:
                not_circ_results.append((circ_intervals, support_mates))
        return circ_results, not_circ_results

    def extend(self, long_interval, other_index_intervals, circ_intervals, support_mates):
        extend_interval = self.interval_extend(long_interval, self.extend_size)
        interval_for_remove = [other_index_interval.other_index_interval for other_index_interval in
                               other_index_intervals]
        interval_for_indexs = [interval_for_index for interval_for_index in self.long_intervals_for_index if
                               interval_for_index not in interval_for_remove]
        intersect_intervals = extend_interval.intersects(interval_for_indexs)
        extend_depth, extend_or_not = self.extend_amplified(long_interval, extend_interval, circ_intervals)
        if intersect_intervals:
            for intersect_interval in intersect_intervals:
                if intersect_interval.clipped == utils.opposite[long_interval.clipped]:
                    if self.similar_depth(long_interval, intersect_interval, circ_intervals):
                        whole_interval = self.combine(long_interval, intersect_interval)
                        self.left_right_depth(whole_interval)
                        if intersect_interval in other_index_intervals:
                            intersect_id = other_index_intervals.index(intersect_interval)
                            yield (
                                None, None, circ_intervals[intersect_id:] + [whole_interval],
                                support_mates[intersect_id:] + [intersect_interval.mate], intersect_id)
                        else:
                            yield (
                                self.shorten(intersect_interval.other_long_interval),
                                other_index_intervals + [intersect_interval],
                                circ_intervals + [whole_interval],
                                support_mates + [intersect_interval.mate], None)
        if extend_or_not:
            if long_interval.clipped == 'left':
                long_interval.end = extend_interval.end
            else:
                long_interval.start = extend_interval.start
            long_interval.extend_depth.append(extend_depth)
            yield (long_interval, other_index_intervals, circ_intervals, support_mates, None)
        else:
            if long_interval.clipped == "right":
                whole_interval = utils.WholeInterval(chrom=long_interval.chrom, start=long_interval.start,
                                                     end=long_interval.end, strand='-',
                                                     depths=long_interval.extend_depth)
            else:
                whole_interval = utils.WholeInterval(chrom=long_interval.chrom, start=long_interval.start,
                                                     end=long_interval.end, strand='+',
                                                     depths=long_interval.extend_depth)
            # avoid influence
            if len(long_interval.extend_depth) > 1:
                long_interval.extend_depth = long_interval.extend_depth[:-1]
            self.left_right_depth(whole_interval)
            yield (None, None, circ_intervals + [whole_interval], support_mates, None)

    def shorten(self, long_interval):
        if long_interval.clipped == 'left':
            long_interval.end = long_interval.start + self.max_insert
        if long_interval.clipped == 'right':
            long_interval.start = long_interval.end - self.max_insert
        return long_interval

    def similar_depth(self, interval1, intersect_interval, previous_intervals):
        extend_depth = interval1.extend_depth
        low, high, low_value, high_value = self.limit(extend_depth, previous_intervals)
        if low <= low_value and low <= intersect_interval.depth < high and high_value <= high and low <= intersect_interval.other_index_interval.depth < high:
            return True
        else:
            return False

    def limit(self, extend_depth, previous_intervals=None):
        sum_depth = [(np.average(extend_depth), len(extend_depth), np.min(extend_depth), np.max(extend_depth))]
        if previous_intervals:
            sum_depth.extend([(interval.depth, len(interval.depths), np.min(interval.depths),
                               np.max(interval.depths)) for interval in previous_intervals])
        low_value = np.min([i[2] for i in sum_depth])
        high_value = np.max([i[3] for i in sum_depth])
        average = np.sum([i[0] * i[1] for i in sum_depth]) / np.sum([i[1] for i in sum_depth])
        low = average - (average / self.depth_average) * self.depth_std
        high = average + (average / self.depth_average) * self.depth_std
        return low, high, low_value, high_value

    def extend_amplified(self, long_interval, extend_interval, previous_intervals):
        extend_part_depth = self.interval_depth(extend_interval)
        if extend_part_depth < self.upper_limit or long_interval.length() > 20000000:
            return extend_part_depth, False
        low, high, low_value, high_value = self.limit(long_interval.extend_depth, previous_intervals)
        if low < low_value and low <= long_interval.depth and low <= extend_part_depth:
            return extend_part_depth, True
        else:
            return extend_part_depth, False

    @staticmethod
    def interval_extend(long_interval, extend_size):
        if long_interval.length() < 10000:
            if long_interval.clipped == 'left':
                return utils.Interval(long_interval.chrom, long_interval.start,
                                      long_interval.end + extend_size)
            else:
                return utils.Interval(long_interval.chrom, long_interval.start - extend_size, long_interval.end)

        else:
            if long_interval.clipped == 'left':
                return utils.Interval(long_interval.chrom, long_interval.end,
                                      long_interval.end + extend_size)
            else:
                return utils.Interval(long_interval.chrom, long_interval.start - extend_size,
                                      long_interval.start)

    @staticmethod
    def combine(current_interval, intersect_interval):
        if current_interval.clipped == "left" and intersect_interval.clipped == 'right':
            interval = utils.WholeInterval(current_interval.chrom, current_interval.start, intersect_interval.end,
                                           strand='+', depths=current_interval.extend_depth)
            return interval
        if current_interval.clipped == "right" and intersect_interval.clipped == 'left':
            interval = utils.WholeInterval(current_interval.chrom, intersect_interval.start, current_interval.end,
                                           strand='-', depths=current_interval.extend_depth)
            return interval

    def interval_depth(self, interval):
        if interval.start <= 0:
            return 0
        else:
            depth = sum([sum(a) for a in self.bam.count_coverage(interval.chrom, interval.start, interval.end)]) / (
                    interval.end - interval.start)
            return depth

    def build_long_intervals(self):
        long_interval1s = [utils.LongInterval(i.interval1, i.clipped1, i) for i in self.split_read_mates]
        long_interval2s = [utils.LongInterval(i.interval2, i.clipped2, i) for i in self.split_read_mates]
        for long_interval in long_interval1s:
            long_interval.depth = long_interval.mate.interval1.depth
            long_interval.extend_depth.append(long_interval.depth)
        for long_interval in long_interval2s:
            long_interval.depth = long_interval.mate.interval2.depth
            long_interval.extend_depth.append(long_interval.depth)
        for long_interval1, long_interval2 in zip(long_interval1s, long_interval2s):
            long_interval1_for_index = copy.copy(long_interval1)
            long_interval2_for_index = copy.copy(long_interval2)

            long_interval1.other_long_interval = long_interval2
            long_interval1.other_index_interval = long_interval2_for_index
            long_interval2_for_index.other_long_interval = long_interval1
            long_interval2_for_index.other_index_interval = long_interval1_for_index

            long_interval2.other_long_interval = long_interval1
            long_interval2.other_index_interval = long_interval1_for_index
            long_interval1_for_index.other_long_interval = long_interval2
            long_interval1_for_index.other_index_interval = long_interval2_for_index
            self.long_intervals_for_index.append(long_interval1_for_index)
            self.long_intervals_for_index.append(long_interval2_for_index)
            self.mates.append((long_interval1, long_interval2, True))
        print(f"finally, {len(self.mates)} mates, {len(self.long_intervals_for_index)} intervals")

    def left_right_depth(self, whole_interval):
        left_interval = utils.Interval(whole_interval.chrom, whole_interval.start - utils.extension,
                                       whole_interval.start)
        right_interval = utils.Interval(whole_interval.chrom, whole_interval.end, whole_interval.end + utils.extension)
        whole_interval.left_depth = self.interval_depth(left_interval)
        whole_interval.right_depth = self.interval_depth(right_interval)

    def filter(self):
        split_read_mates = self.split_read_mates
        print(f"before process, {len(split_read_mates)} sr mates")
        if self.cutoff:
            split_read_mates = [split_read_mate for split_read_mate in split_read_mates if
                                len(split_read_mate.support_split_reads) >= self.split_read_cutoff]
            split_read_mates = [split_read_mate for split_read_mate in split_read_mates if
                                len(split_read_mate.support_discordant_reads) >= self.discordant_read_cutoff]
        print(f"filter, {len(split_read_mates)} sr mates")
        removes = []
        for split_read_mate in split_read_mates:
            depth1 = split_read_mate.interval1.depth
            depth2 = split_read_mate.interval2.depth
            low, high, low_value, high_value = self.limit([depth1, depth2])
            if depth1 > self.upper_limit and depth2 > self.upper_limit and low < depth1 < high and low < depth2 < high:
                pass
            else:
                removes.append(split_read_mate)
        split_read_mates = [split_read_mate for split_read_mate in split_read_mates if
                            split_read_mate not in removes]
        print(f"filter, {len(split_read_mates)} sr mates")
        split_read_mates = utils.process_all_sr_mates(split_read_mates)
        print(f"after process, {len(split_read_mates)} sr mates")
        self.split_read_mates = split_read_mates
