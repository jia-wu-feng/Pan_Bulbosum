# Pan_Bulbosum

Code for Article "A haplotype-resolved pangenome of the barley wild relative Hordeum bulbosum"

## Split by chromosome

````shell

for i in chr1H chr2H chr3H chr4H chr5H chr6H chr7H;
do
echo $i
name=`basename $i`
echo $name
mkdir $name
cd $name

date


for n in `cat genome_list.txt`;do echo $n; cname=`basename $n .fasta` ; #/opt/Bio/seqkit/0.9.1/bin/seqkit grep -f 1 $n > $cname.fa ; done

for m in `cat genome_list.txt`;do echo $m;qname=`basename -s .fasta $m`; echo $qname; #sed -i '1c >'$qname'.'$name'' $qname.fa ; done

for m in `ls /filer-dg/agruppen/dg7/fengj/minigraph_cactus/HB_minigraph_cactus/$name/*fa`; do qname=`basename -s .fa $m`; echo -e $qname"\t"$m >> $name.txt; done

cd /filer-dg/agruppen/dg7/fengj/minigraph_cactus/HB_minigraph_cactus 

done

````

## Execute Minigraph-Cactus Pangenome Pipeline

````shell

sfor i in chr1H chr2H chr3H chr4H chr5H chr6H chr7H;
do
echo $i
name=`basename $i`
echo $name
cd $name

echo "#!/bin/bash

date

mkdir cactus-scratch

singularity exec -B /filer-dg/agruppen --pwd \$PWD /filer-dg/agruppen/dg7/fengj/minigraph_cactus/HB_minigraph_cactus/cactus_v2.6.13.sif cactus-pangenome ./js ./$name.txt --outDir ./$name --outName hb.$name --reference FB190113Hap1 --filter 4 --giraffe clip filter --viz --odgi --chrom-vg clip filter --chrom-og --gbz clip filter full --gfa clip full --vcf --vcfReference FB190113Hap1 --logFile ./hb.log --workDir ./cactus-scratch --consCores 40 --mgMemory 500Gi

date " > cactus.zsh

sbatch cactus.zsh

cd /filer-dg/agruppen/dg7/fengj/minigraph_cactus/HB_minigraph_cactus 

done

````
