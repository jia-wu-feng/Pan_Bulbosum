#!/bin/zsh

projectdir="/filer-dg/agruppen/dg2/fengj/hbulbosum_genome/genome_assembly/A40"

reads=$projectdir/hifi_reads

prefix=${projectdir:r}/A40_hifiasm

echo $prefix

find $reads | grep 'fastq.gz$' | xargs /filer-dg/agruppen/dg7/fengj/software/hifiasm-0.13/hifiasm -t 50 -o $prefix  2> ${prefix}_asm.err


