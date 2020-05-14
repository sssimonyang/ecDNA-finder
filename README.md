### prepare
To run the software simply, you only need to give it the path to coordinate-sorted bam file and queryname-sorted bam file, and a directory name.

Use follow command to transform a coordinate-sorted bam file to queryname-sorted.
```shell script
samtools sort -n -@ {threads_num} -o {queryname-sorted} {coordinate-sorted}
```
-------------------------
### run
run simply
```shell script
python main.py -c {coordinate-sorted} -q {queryname-sorted} -d {dirname}
```

### help guide
```
usage: main.py [-h] -c COORDINATE -q QUERYNAME -d DIRNAME

ecDNA-finder

optional arguments:
  -h, --help            show this help message and exit
  -c COORDINATE, --coordinate COORDINATE
                        bam file sorted by coordinate
  -q QUERYNAME, --queryname QUERYNAME
                        bam file sorted by queryname
  -d DIRNAME, --dirname DIRNAME
                        result directory
```

### pbs_file
pbs_file recommended for lab use
```shell script
#!/bin/sh
#PBS -N PBS_ecDNA
#PBS -l nodes=1:ppn=5
#PBS -l walltime=16:00:00
#PBS -S /bin/bash
#PBS -q normal_3
#PBS -o /public/home/zhangjing1/yangjk/ecDNA/result/result_out
#PBS -e /public/home/zhangjing1/yangjk/ecDNA/result/result_error

start='------------START------------'$(date "+%Y %h %d %H:%M:%S")'------------START------------'
echo $start
echo $start >&2
dirname=SRR8236745
cd /public/home/zhangjing1/yangjk/ecDNA/result/
coordinate=/public/home/zhangjing1/yangjk/data/bam/${dirname}/sorted_coordinate.bam
queryname=/public/home/zhangjing1/yangjk/data/bam/${dirname}/sorted_query_name.bam
/public/home/liuxs/anaconda3/envs/ecDNA/bin/python /public/home/zhangjing1/yangjk/ecDNA/code/main.py -c ${coordinate} -q ${queryname} -d ${dirname}
end='-------------END-------------'$(date "+%Y %h %d %H:%M:%S")'-------------END-------------'
echo $end
echo $end >&2

```

