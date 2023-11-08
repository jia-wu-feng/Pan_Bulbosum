#!/bin/zsh

bcftools='/opt/Bio/bcftools/1.12/bin/bcftools'
bgzip='/opt/Bio/bcftools/1.12/bin/bgzip'
tabix='/opt/Bio/bcftools/1.12/bin/tabix'

read <<< $* ref bamlist prefix region part

base=${prefix}_${part}_${region}

exec 2> $base.err 

echo $PWD 1>&2
echo $HOST 1>&2
echo $SLURM_JOB_ID 1>&2
echo $* 1>&2
date 1>&2
$bcftools mpileup --threads 2 -d 30000000 -q 20 -a DP,AD -b $bamlist -r $region -f $ref \
 | $bcftools call --threads 2 -mv - \
 | $bgzip -@ 2 > $base.vcf.gz 

echo $pipestatus | tee ${base}_pipestatus.txt  | tr ' ' '\n' | grep -q '^[^0$]' || \
 $tabix -Cp vcf $base.vcf.gz 
date 1>&2
