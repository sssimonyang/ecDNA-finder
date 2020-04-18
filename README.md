### prepare
Make a directory
The directory needs a coordinate-sorted bam file named `sorted_coordinate.bam` and a queryname=sorted bam file named `sorted_query_name.bam`
Use follow commands to tranform a coordinate-sorted bam file to queryname-sorted
```shell script
samtools sort -@ {threads_num} -o sorted_coordinate.bam sorted_coordinate.bam
```

### run
change working directory and run
```shell script
/public/home/liuxs/anaconda3/envs/ecDNA/bin/python /public/home/zhangjing1/yangjk/ecDNA/code/main.py
```


### pbs_file
```shell script
#!/bin/sh
#PBS -N PBS_main
#PBS -l nodes=1:ppn=5
#PBS -l walltime=24:00:00
#PBS -S /bin/bash
#PBS -q normal_3
#PBS -o /public/home/zhangjing1/yangjk/pbs_main_out
#PBS -e /public/home/zhangjing1/yangjk/pbs_main_error
cd /public/home/zhangjing1/yangjk/data
sraid=SRR8236746
start='------------START------------'$(date "+%Y %h %d %H:%M:%S")'------------START------------'
echo $start
echo $start >&2
samtools sort -n -@ 24 -o sorted_query_name.bam sorted_coordinate.bam
/public/home/liuxs/anaconda3/envs/ecDNA/bin/python /public/home/zhangjing1/yangjk/ecDNA/code/main.py
end='-------------END-------------'$(date "+%Y %h %d %H:%M:%S")'-------------END-------------'
echo $end
echo $end >&2
```

