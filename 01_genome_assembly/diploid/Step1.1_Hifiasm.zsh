#!/bin/zsh

projectdir="/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/FB20_005_1"

reads=$projectdir/hifi_reads

prefix=${projectdir:r}/FB20_005_1_hifiasm

echo $prefix

/home/fengj/tools/hifiasm-0.13/hifiasm -v

find $reads | grep 'fastq.gz$' | xargs /home/fengj/tools/hifiasm-0.13/hifiasm -t 50 -o $prefix  2> ${prefix}_asm.err
