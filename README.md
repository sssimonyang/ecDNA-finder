This software is a modification and improvement of [CircleMap](https://github.com/iprada/Circle-Map) which is under MIT license to detect ecDNA from WGS bam file.

extract.py utils.py mates.py is an modification of CircleMap/circlemap/extract_circle_SV_reads.py CircleMap/circlemap/utils.py CircleMap/circlemap/repeats.py respectively.

Thanks to @iprada for developing circelmap and facilitate my development of ecDNA-finder.


### Prepare
To run the software simply, you only need to prepare coordinate-sorted bam file and queryname-sorted bam file.

Use follow command to transform a coordinate-sorted bam file to queryname-sorted.
```shell script
samtools sort -n -@ {threads_num} -o {queryname-sorted} {coordinate-sorted}
```

### Package need
- numpy
- pandas
- pysam

> All can be easy installed by pip or conda

### Run
Run simply
```shell script
python main.py -coord {coordinate-sorted} -query {queryname-sorted} -dir {dirname}
```

### Other important params you may need
-cutoff : the default value of cutoff is 0, which means use the round(depth_average / 20) to cutoff peak and split read mate.
The depth_average is calculated automatically. Certainly, the mininum allowable value is 1. YOU can provide this value to change the result.
BUT the amount of time used varies accordingly.

### Results
Results is placed in a new directory named {dirname} in the current directory.
The circ_results.tsv is a tab-delimiter file.

### Help guide
```
usage: main.py [-h] -coord COORDINATE -query QUERYNAME -dir DIRNAME
                   [-cutoff CUTOFF]

ecDNA-finder

optional arguments:
  -h, --help            show this help message and exit
  -coord COORDINATE, --coordinate COORDINATE
                        bam file sorted by coordinate
  -query QUERYNAME, --queryname QUERYNAME
                        bam file sorted by queryname
  -dir DIRNAME, --dirname DIRNAME
                        result directory
  -cutoff CUTOFF, --cutoff CUTOFF
                        seed interval cutoff and support read cutoff
```

### Pbs_file
A pbs_file recommended for lab use
```shell script
#!/bin/sh
#PBS -N PBS_ecDNA
#PBS -l nodes=1:ppn=5
#PBS -l walltime=16:00:00
#PBS -S /bin/bash
#PBS -q normal_3
#PBS -o /public/home/zhangjing1/yangjk/ecDNA/result/ecDNA_out
#PBS -e /public/home/zhangjing1/yangjk/ecDNA/result/ecDNA_out

start='------------START------------'$(date "+%Y %h %d %H:%M:%S")'------------START------------'
echo $start
echo $start >&2
dirname=SRR8236745
cd /public/home/zhangjing1/yangjk/ecDNA/result/
coordinate=/public/home/zhangjing1/yangjk/data/bam/${dirname}/sorted_coordinate.bam
queryname=/public/home/zhangjing1/yangjk/data/bam/${dirname}/sorted_query_name.bam
/public/home/liuxs/anaconda3/envs/ecDNA/bin/python /public/home/zhangjing1/yangjk/ecDNA/code/main.py -coord ${coordinate} -query ${queryname} -dir ${dirname}
end='-------------END-------------'$(date "+%Y %h %d %H:%M:%S")'-------------END-------------'
echo $end
echo $end >&2

```

