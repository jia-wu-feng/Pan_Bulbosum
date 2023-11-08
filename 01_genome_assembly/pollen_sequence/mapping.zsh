#!/bin/zsh

adapter=$1
shift
minlen=$1
shift
threads=$1
shift
novosort_threads=$1
shift
mem=$1
shift
ref=$1
shift
tmp=$1
shift
i=$1

minimap='/opt/Bio/minimap2/2.17/bin/minimap2'
samtools='/opt/Bio/samtools/1.10/bin/samtools'
novosort="/opt/Bio/novocraft/V3.06.05/bin/novosort"
cutadapt='/opt/Bio/cutadapt/3.3/bin/cutadapt'

indexsize='50G'
batchmem='5G'

base=$i/$i:t
bam=${base}.bam

minimaperr=${base}_minimap.err
samtoolserr=${base}_samtools.err
sorterr=${base}_novosort.err
cutadapterr=${base}_cutadapt.err

b=$base:t
rgentry="@RG\tID:$b\tPL:ILLUMINA\tPU:$b\tSM:$b"

$cutadapt -j $threads --interleaved -a $adapter -A $adapter -m $minlen \
 <(find $i | egrep '_R1_0' | grep 'f.*q.gz$' | sort | xargs gzip -cd) \
 <(find $i | egrep '_R2_0' | grep 'f.*q.gz$' | sort | xargs gzip -cd) \
 2> $cutadapterr \
 | $minimap -ax sr -R $rgentry -t $threads -2 -I $indexsize -K$batchmem -a $ref /dev/stdin 2> $minimaperr \
 | $samtools view -Su /dev/stdin 2> $samtoolserr \
 | $novosort -c $novosort_threads -t $tmp -m $mem -i --keepTags --md -o $bam /dev/stdin 2> $sorterr

echo $pipestatus > ${base}_pipestatus.txt 
