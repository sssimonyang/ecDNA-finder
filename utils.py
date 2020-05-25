#!/usr/bin/env python 
# -*- coding: utf-8 -*-
"""
@author: SSSimon Yang
@contact: yangjingkang@126.com
@file: utils.py
@time: 2020/3/11 20:03
@desc:
"""
import itertools
import itertools as it
from collections import namedtuple

import numpy as np
import pandas as pd
import pysam as ps

Read = namedtuple("Read",
                  ['query_name', 'reference_name', 'reference_start', 'reference_end', 'strand', 'query_sequence'])

ReadMate = namedtuple("ReadMate",
                      ['type', 'read1', 'read2'])
opposite = {'right': 'left', 'left': 'right'}
# cutoff
mapq_cutoff = 40
interval_p_cutoff = 0.5

# number
sample_size = 10000
bed_size = 10000
std = 4

# bp
extend_size = 10000
cluster_distance = 200
tlen_min = 5000

# calculate automaticly or infer from the given bam
cores = 4
extension = 900
insert_mean = 500
insert_std = 100
depth_average = 60
depth_std = 20

cutoff = 0
chrom_names = []
chrom_lengths = []


def colinear(read, read_mate, supplementary):
    if clipped_direction(read) in ('right', 'left') and clipped_direction(supplementary) in ('right', 'left'):
        if read.is_reverse == supplementary.is_reverse and \
                clipped_direction(read) == opposite[clipped_direction(supplementary)] or \
                read.is_reverse != supplementary.is_reverse and \
                clipped_direction(read) == clipped_direction(supplementary):
            if is_close(read, read_mate) and not is_close(read, supplementary):
                # if (clipped_direction(read) == 'right' and read_mate.reference_start < read.reference_start) or \
                #         (clipped_direction(read) == 'left' and read_mate.reference_start > read.reference_start):
                return True
            if not is_close(read, read_mate) and is_close(read_mate, supplementary):
                # if (clipped_direction(supplementary) == 'right' and read_mate.reference_start <
                # supplementary.reference_start) or \ (clipped_direction(supplementary) == 'left' and
                # read_mate.reference_start > supplementary.reference_start):
                return True
    return False


def is_close(read1, read2):
    if read1.reference_name != read2.reference_name:
        return False
    else:
        if abs(read1.reference_start - read2.reference_start) >= tlen_min:
            return False
        else:
            return True


def is_clipped(read):
    cigar = [''.join(g) for _, g in it.groupby(read.cigarstring, str.isalpha)]
    if "H" in cigar or 'S' in cigar:
        return True
    else:
        return False


def clipped_type(read):
    cigar = [''.join(g) for _, g in it.groupby(read.cigarstring, str.isalpha)]
    if "H" in cigar:
        return 'hard'
    if "S" in cigar:
        return 'soft'


def clipped_direction(read):
    cigar = [''.join(g) for _, g in it.groupby(read.cigarstring, str.isalpha)]
    if cigar[1] in ['S', 'H'] and cigar[-1] in ['S', 'H'] and len(cigar) >= 6:
        return 'both'
    if cigar[1] in ['S', 'H'] and len(cigar) >= 4:
        return 'left'
    if cigar[-1] in ['S', 'H'] and len(cigar) >= 4:
        return 'right'


class SA:
    def __init__(self, info):
        chrom, start, strand, cigarstring, mapping_quality, edit_distance = [
            i.strip() for i in info.split(";")[0].split(',')]
        self.chrom = chrom
        self.start = int(start)
        self.strand = strand
        self.cigarstring = cigarstring
        self.mapping_quality = int(mapping_quality)
        self.edit_distance = int(edit_distance)
        self.is_reverse = True if self.strand == "-" else False
        self.cigar_list = [''.join(i) for _, i in it.groupby(
            self.cigarstring, str.isalpha)]
        self.end = self.start + self.genome_alignment()

    def genome_alignment(self):
        """Function that gets as input the SA tag CIGAR and returns the length of the alignment interval in the
        genome it will look at the number of matches and deletions in the CIGAR, as they are the elements that will
        explain the genome alignment
        """
        aligned = 0

        # do it for the matches
        match_index = [x for x in range(
            len(self.cigar_list)) if self.cigar_list[x] == 'M']
        for index in match_index:
            aligned += int(self.cigar_list[index - 1])

        match_index = [x for x in range(
            len(self.cigar_list)) if self.cigar_list[x] == 'D']
        for index in match_index:
            aligned += int(self.cigar_list[index - 1])

        return aligned

    def __getattr__(self, item):
        if item == 'reference_name':
            return self.chrom
        if item == 'reference_start':
            return self.start


def show_mates(reads, sr):
    if sr:
        return [ReadMate('SR', Read(read.query_name, read.reference_name, read.reference_start, read.reference_end,
                                    "-" if read.is_reverse else '+',
                                    read.query_sequence),
                         Read(read.query_name, mate.chrom, mate.start, mate.end,
                              "-" if mate.is_reverse else '+',
                              read.query_sequence)) for read, mate in reads]
    else:
        return [ReadMate('DR', Read(read.query_name, read.reference_name, read.reference_start, read.reference_end,
                                    "-" if read.is_reverse else '+',
                                    read.query_sequence),
                         Read(mate.query_name, mate.reference_name, mate.reference_start, mate.reference_end,
                              "-" if mate.is_reverse else '+',
                              mate.query_sequence)) for read, mate in reads]


def phred_to_prob(values):
    """Function that takes as input a numpy array with phred base quality scores and returns an array with base probabi-
    lity scores"""
    return 10 ** ((values * -1) / 10)


class SR_Mate:
    def __init__(self, interval1, clipped1, interval2, clipped2):
        self.interval1 = interval1
        self.interval2 = interval2
        self.clipped1 = clipped1
        self.clipped2 = clipped2
        self.support_discordant_reads = []
        self.support_split_reads = []
        self.from_read = None
        self.support_discordant_mates = []

    def __repr__(self):
        return ' '.join(['SR_Mate', str(self.interval1), self.clipped1, 'and', str(self.interval2), self.clipped2])

    __str__ = __repr__

    # deprecated
    def out(self):
        return [str(self.interval1), self.clipped1, str(self.interval2), self.clipped2,
                ' '.join([read.query_name for read in self.support_split_reads]),
                ' '.join([read.query_name for read in self.support_discordant_reads]),
                ' '.join([str(discordant_mate) for discordant_mate in self.support_discordant_mates])]

    def get_points(self):
        points = []
        if self.clipped1 == 'left':
            points.append(Point(self.interval1.chrom, self.interval1.start, 'left'))
        else:
            points.append(Point(self.interval1.chrom, self.interval1.end, 'right'))
        if self.clipped2 == 'left':
            points.append(Point(self.interval2.chrom, self.interval2.start, 'left'))
        else:
            points.append(Point(self.interval2.chrom, self.interval2.end, 'right'))
        return points


class Discordant_Mate:
    def __init__(self, interval1, interval2):
        self.interval1 = interval1
        self.interval2 = interval2
        self.support_discordant_reads = []

    def __repr__(self):
        return ' '.join(['Discordant_Mate', str(self.interval1), 'and', str(self.interval2)])

    __str__ = __repr__

    # deprecated
    def out(self):
        return [str(self.interval1), str(self.interval2),
                ' '.join([read.query_name for read in self.support_discordant_reads])]


def process_sr_mate(in_sr_mates):
    results = []
    left_sr_mates = [sr_mate for sr_mate in in_sr_mates if sr_mate.clipped1 == 'left']
    right_sr_mates = [sr_mate for sr_mate in in_sr_mates if sr_mate.clipped1 == 'right']
    for sr_mates in [left_sr_mates, right_sr_mates]:
        if not sr_mates:
            continue
        interval1 = sr_mates[0].interval1
        clipped1 = sr_mates[0].clipped1
        interval1.start = min([i.start for i in [j.interval1 for j in sr_mates]])
        interval1.end = max([i.end for i in [j.interval1 for j in sr_mates]])

        left = [i.interval2 for i in sr_mates if i.clipped2 == 'left']
        right = [i.interval2 for i in sr_mates if i.clipped2 == 'right']
        if left:
            interval2 = left[0]
            interval2.start = min([i.start for i in left])
            interval2.end = max([i.end for i in left])
            result = SR_Mate(interval1, clipped1, interval2, 'left')
            result.support_split_reads = [i.from_read for i in sr_mates if i.clipped2 == 'left']
            results.append(result)
        if right:
            interval2 = right[0]
            interval2.start = min([i.start for i in right])
            interval2.end = max([i.end for i in right])
            result = SR_Mate(interval1, clipped1, interval2, 'right')
            result.support_split_reads = [i.from_read for i in sr_mates if i.clipped2 == 'right']
            results.append(result)
    return results


class Interval:
    def __init__(self, chrom, start, end, strand=''):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        assert strand in ['+', '-', '']
        self.strand = strand
        self.depth = None

    def contain(self, interval):
        if self.chrom == interval.chrom and self.start <= interval.start and self.end >= interval.end:
            return True
        else:
            return False

    def length(self):
        return self.end - self.start

    def __repr__(self):
        return self.__class__.__name__ + ' ' + ','.join([self.chrom, str(self.start), str(self.end)])

    __str__ = __repr__

    def intersect(self, interval):
        if self.chrom != interval.chrom:
            return False
        if self.end < interval.start:
            return False
        if self.start > interval.end:
            return False
        return True

    def intersects(self, intervals):
        intersect_interval = []
        for interval in intervals:
            if self.intersect(interval):
                intersect_interval.append(interval)
        return intersect_interval


def to_interval(interval):
    if isinstance(interval, pd.Series):
        return Interval(chrom=interval.chrom, start=interval.start, end=interval.end)
    if isinstance(interval, ps.libcalignedsegment.AlignedSegment):
        out = Interval(chrom=interval.reference_name, start=interval.reference_start, end=interval.reference_end,
                       strand='-' if interval.is_reverse else '+')
        return out
    if isinstance(interval, list):
        return Interval(chrom=interval[0], start=interval[1], end=interval[2])


def get_mate_intervals(bam, interval, mapq_cutoff=40):
    candidate_mates = []
    for read in bam.fetch(interval.chrom, interval.start, interval.end, multiple_iterators=True):

        if read.mapping_quality >= mapq_cutoff:

            if clipped_type(read) != 'hard' and (read.tlen == 0 or read.tlen > tlen_min):
                read_mate = bam.mate(read)
                mate_interval = [read_mate.reference_name, read_mate.reference_start,
                                 read_mate.reference_end, "DR",
                                 'Both',
                                 str(1 - phred_to_prob(read_mate.mapping_quality))]
                candidate_mates.append(mate_interval)

            if is_clipped(read) and read.has_tag('SA'):
                supplementary = SA(read.get_tag('SA'))

                if supplementary.mapping_quality >= mapq_cutoff:
                    if (read.is_reverse == supplementary.is_reverse and clipped_direction(read) == opposite[
                        clipped_direction(supplementary)]) or \
                            (read.is_reverse != supplementary.is_reverse and clipped_direction(
                                read) == clipped_direction(supplementary)):
                        alignment_length = supplementary.genome_alignment()

                        if clipped_direction(supplementary) == 'left':
                            mate_interval = [supplementary.chrom, supplementary.start - alignment_length,
                                             supplementary.end, "SA", opposite[clipped_direction(supplementary)],
                                             str(1 - phred_to_prob(supplementary.mapping_quality))]
                        else:
                            mate_interval = [supplementary.chrom, supplementary.start,
                                             supplementary.end + alignment_length, "SA",
                                             opposite[clipped_direction(supplementary)],
                                             str(1 - phred_to_prob(supplementary.mapping_quality))]
                        candidate_mates.append(mate_interval)

    return candidate_mates


def get_realignment_intervals(intervals, interval_p_cutoff=0.5):
    labels = ['chrom', 'start', 'end', 'read_type', 'orientation', 'probability']
    candidate_mates_dataframe = pd.DataFrame.from_records(intervals, columns=labels)
    candidate_mates_dataframe = candidate_mates_dataframe.sort_values(by=['chrom', 'start', 'end'],
                                                                      ascending=[True, True, True])
    candidate_mates_dataframe['probability'] = candidate_mates_dataframe.probability.astype(float)
    candidate_mates = candidate_mates_dataframe.groupby(
        ((candidate_mates_dataframe.end.shift() - candidate_mates_dataframe.start).lt(-cluster_distance) | (
                candidate_mates_dataframe.chrom.shift() != candidate_mates_dataframe.chrom)
         ).cumsum()).agg(
        {'chrom': 'first', 'start': 'first', 'end': 'last', 'probability': 'sum'})
    candidate_mates['probability'] = candidate_mates['probability'] / candidate_mates.probability.sum()
    if interval_p_cutoff == 0:
        candidate_mates = candidate_mates.loc[
            candidate_mates['probability'] == candidate_mates['probability'].max()]
    else:
        candidate_mates = candidate_mates.loc[
            candidate_mates['probability'] >= interval_p_cutoff]
    return candidate_mates.sort_values(by=['probability'], ascending=[False])


def assign_discordant_reads(split_read_mates, discordant_reads):
    for split_read_mate in split_read_mates:
        for read1, read2 in discordant_reads:
            read1_interval = to_interval(read1)
            read2_interval = to_interval(read2)
            if split_read_mate.interval1.contain(read1_interval) and split_read_mate.interval2.contain(read2_interval):
                split_read_mate.support_discordant_reads.append(show_mates([(read1, read2)], sr=False)[0])
            elif split_read_mate.interval1.contain(read2_interval) and split_read_mate.interval2.contain(
                    read1_interval):
                split_read_mate.support_discordant_reads.append(show_mates([(read1, read2)], sr=False)[0])
    return split_read_mates


def assign_discordant_mates(split_read_mates, discordant_mates):
    for split_read_mate in split_read_mates:
        for discordant_mate in discordant_mates:
            read1 = discordant_mate.interval1
            read2 = discordant_mate.interval2
            if (split_read_mate.interval1.contain(read1) and split_read_mate.interval2.contain(read2)) or \
                    (split_read_mate.interval2.contain(read1) and split_read_mate.interval1.contain(read2)):
                for read in discordant_mate.support_discordant_reads:
                    split_read_mate.support_discordant_reads.append(read)
                split_read_mate.support_discordant_mates.append(discordant_mate)
    remove_discordant_mates = [discordant_mate for split_read_mate in split_read_mates for discordant_mate in
                               split_read_mate.support_discordant_mates]
    discordant_mates = [discordant_mate for discordant_mate in discordant_mates
                        if discordant_mate not in remove_discordant_mates]
    return split_read_mates, discordant_mates


def get_longest_soft_clipped_bases(reads):
    results = []
    for read in reads:
        cigar = [''.join(g) for _, g in it.groupby(read.cigarstring, str.isalpha)]
        if cigar[1] in ['S', 'H'] and cigar[-1] in ['S', 'H'] and len(cigar) >= 6:
            if int(cigar[0]) > int(cigar[-2]):
                results.append(('left', int(cigar[0]), read))
            else:
                results.append(('right', int(cigar[-2]), read))
        if cigar[1] in ['S', 'H'] and len(cigar) >= 4:
            results.append(('left', int(cigar[0]), read))
        if cigar[-1] in ['S', 'H'] and len(cigar) >= 4:
            results.append(('right', int(cigar[-2]), read))

    results.sort(key=lambda x: x[1], reverse=True)
    return results[0]


class LongInterval(Interval):
    def __init__(self, interval, clipped, mate):
        self.interval = interval
        super().__init__(interval.chrom, interval.start, interval.end)
        self.clipped = clipped
        self.mate = mate
        self.raw_depth = 0.0
        self.other_long_interval = None
        self.other_index_interval = None
        self.extend_depth = []


class WholeInterval(Interval):
    def __init__(self, chrom, start, end, strand, depths):
        super().__init__(chrom, start, end, strand)
        self.depths = depths
        self.left_depth = None
        self.right_depth = None
        self.depth = np.average(self.depths)


def circ_result_out(circ_results):
    out = []
    if not len(circ_results):
        return
    for i, circ in enumerate(circ_results):
        for interval, mate in zip(circ[0], circ[1]):
            out.append(
                [interval.chrom, interval.start, interval.end, interval.strand, interval.length(),
                 interval.depth, interval.left_depth, interval.right_depth,
                 len(mate.support_split_reads),
                 len(mate.support_discordant_reads), i + 1])

    out = pd.DataFrame(out, columns=['chrom', 'start', 'end', 'strand', 'length', 'average_depth', 'left_depth',
                                     'right_depth', 'support_split_reads', 'support_discordant_reads', 'circ_id'])
    group = out.groupby('circ_id').aggregate(
        {'circ_id': 'count', 'average_depth': 'mean', 'chrom': 'unique', 'length': 'sum'})
    group.sort_index()
    group.drop_duplicates(subset=['circ_id', 'length'],
                          keep='first', inplace=True)
    out = out[out.circ_id.isin(group.index)]
    replace = group.reset_index(drop=True).index
    to_replace = group.index
    out['circ_id'] = out['circ_id'].replace(to_replace, replace + 1)
    print(f'Write success circ {out.circ_id.nunique()}')
    out.to_csv('circ_results.tsv', sep='\t', index=False)


def not_circ_result_out(not_circ_results):
    out = []
    if not len(not_circ_results):
        return
    for i, not_circ in enumerate(not_circ_results):
        not_circ_interval = not_circ[0][-1]
        circ_interval = not_circ[0][:-1]
        circ_support = not_circ[1]
        assert len(circ_interval) == len(circ_support)
        for interval, mate in zip(circ_interval, circ_support):
            out.append(
                [interval.chrom, interval.start, interval.end, interval.strand, interval.length(),
                 interval.depth, interval.left_depth, interval.right_depth,
                 len(mate.support_split_reads),
                 len(mate.support_discordant_reads), i + 1])
        interval = not_circ_interval
        out.append([interval.chrom, interval.start, interval.end, interval.strand, interval.length(),
                    interval.depth, interval.left_depth, interval.right_depth,
                    0, 0, i + 1])

    out = pd.DataFrame(out, columns=['chrom', 'start', 'end', 'strand', 'length', 'average_depth', 'left_depth',
                                     'right_depth', 'support_split_reads', 'support_discordant_reads', 'circ_id'])
    out.to_csv('not_circ_results.tsv', sep='\t', index=False)


class Point:
    def __init__(self, chrom, point, clipped):
        self.chrom = chrom
        self.point = point
        self.clipped = clipped

    def similar(self, point):
        if self.chrom == point.chrom and abs(self.point - point.point) <= 10 and self.clipped == point.clipped:
            return True
        else:
            return False

    def __repr__(self):
        return self.__class__.__name__ + ' ,'.join([self.chrom, str(self.point), str(self.clipped)])

    __str__ = __repr__


def process_all_sr_mates(split_read_mates):
    removes = []
    for sr_mate1, sr_mate2 in itertools.combinations(split_read_mates, 2):
        if (
                sr_mate1.interval1.chrom == sr_mate2.interval1.chrom and sr_mate1.interval2.chrom == sr_mate2.interval2.chrom) or (
                sr_mate1.interval2.chrom == sr_mate2.interval1.chrom and sr_mate1.interval1.chrom == sr_mate2.interval2.chrom):
            if (
                    sr_mate1.clipped1 == sr_mate2.clipped1 and sr_mate1.clipped2 == sr_mate2.clipped2) or (
                    sr_mate1.clipped2 == sr_mate2.clipped1 and sr_mate1.clipped1 == sr_mate2.clipped2):
                points1 = sr_mate1.get_points()
                points2 = sr_mate2.get_points()
                if (points1[0].similar(points2[1]) and points1[1].similar(points2[0])) or (
                        points1[0].similar(points2[0]) and points1[1].similar(points2[1])):
                    removes.append(sr_mate2)
    split_read_mates = [split_read_mate for split_read_mate in split_read_mates
                        if split_read_mate not in removes]
    return split_read_mates


# on 7125 split_read_mates, cost time decrease from 81s to 2s
def process_all_sr_mates(split_read_mates):
    intervals = [
        [[i, split_read_mate.interval1.chrom, split_read_mate.interval1.start, split_read_mate.interval1.end],
         [i, split_read_mate.interval2.chrom, split_read_mate.interval2.start, split_read_mate.interval2.end]]
        for i, split_read_mate in enumerate(split_read_mates)]
    intervals = [j for i in intervals for j in i]
    intervals = pd.DataFrame.from_records(intervals, columns=['sr_id', 'chrom', 'start', 'end'])
    intervals = intervals.sort_values(by=['chrom', 'start', 'end'],
                                      ascending=[True, True, True])
    intervals = intervals.reset_index(drop=True)
    indexs = intervals[(intervals.end.shift() - intervals.start).ge(0) & (
            intervals.chrom.shift() == intervals.chrom)].index
    pairs = [(intervals.loc[index].sr_id, intervals.loc[index - 1].sr_id) for index in indexs]
    pairs = [(sr_id1, sr_id2) if sr_id1 < sr_id2 else (sr_id2, sr_id1) for sr_id1, sr_id2 in pairs]
    pairs = list(set(pairs))
    removes = [split_read_mates[sr_id2] for sr_id1, sr_id2 in pairs if
               similar_check(split_read_mates[sr_id1], split_read_mates[sr_id2])]
    split_read_mates = [split_read_mate for split_read_mate in split_read_mates
                        if split_read_mate not in removes]
    return split_read_mates


def similar_check(split_read_mate1, split_read_mate2):
    points1 = split_read_mate1.get_points()
    points2 = split_read_mate2.get_points()
    if (points1[0].similar(points2[1]) and points1[1].similar(points2[0])) or (
            points1[0].similar(points2[0]) and points1[1].similar(points2[1])):
        return True
    else:
        return False
