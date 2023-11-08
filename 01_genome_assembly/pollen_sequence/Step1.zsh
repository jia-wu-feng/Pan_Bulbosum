#working directory: /filer-dg/agruppen/dg6/mascher/DG/hb_pollenseq_210503

#Set up directory strucuture and copy scripts

d='/hsm/novaseq/CSF/NOVASEQ/2021/210430_A00550_0172_AHYTYNDRXX/houben'

find $d -maxdepth 1 -mindepth 1 -type d | while read i; do
 echo $i:t | cut -d - -f 1 | sed 's/^/Sample_/' | read dd
 mkdir -p $dd
 find $i -type f | xargs ln -st $dd 
done

### use stimulating reads for haplotype2 and haplotype1 
wgsim -e 0 -1 150 -2 150 -r 0 -R 0 -X 0 -N 350000000 FB19_011_3_pseudomolecules_v6_hap1.fasta haplotype1.read1.fq haplotype1.read2.fq

wgsim -e 0 -1 150 -2 150 -r 0 -R 0 -X 0 -N 350000000 FB19_011_3_pseudomolecules_v6_hap2.fasta haplotype2.read1.fq haplotype2.read2.fq

#run the mapping for all sample and stimulating reads
cd /filer-dg/agruppen/DG/mascher/hb_pollenseq_210503
adapter='AGATCGGAAGAGC'
minlen=30
threads=40
novosort_threads=16
mem='20G'
tmp=$TMPDIR
ref='/filer-dg/agruppen/dg2/fengj/hbulbosum_genome/SNP_shortreads/ref_genome/FB19_011_3_pseudomolecules_v6_hap1.fasta'

find -mindepth 1 -type d | parallel -j 100 srun --auks=yes -c 40 --mem '100G' \
 ./mapping.zsh $adapter $minlen $threads $novosort_threads $mem $ref $tmp '{}' &

#get csi instead of bai indices
find | grep 'bai' | grep Sample_24 | xargs rm -f

find | grep 'bam$' | grep Sample_24 | parallel -j 100 srun --auks=yes -c 1 --mem '1G' \
 /opt/Bio/samtools/1.9/bin/samtools index -c '{}'


adapter='AGATCGGAAGAGC'
minlen=30
threads=4
novosort_threads=4
mem='20G'
tmp=$TMPDIR
ref='/filer-dg/agruppen/dg2/fengj/hbulbosum_genome/SNP_shortreads/ref_genome/FB19_011_3_pseudomolecules_v6_hap1.fasta'
srun --auks=yes -c 10 -p gpu --mem '100G' ./mapping.zsh $adapter $minlen $threads $novosort_threads $mem $ref $tmp 


#get BAM list pollen sequence + stimulating reads for haplotype2 and haplotype1
find -L $PWD | grep 'bam$' > 210504_bam_list.txt

#split reference into 5 Mb bins
ref='/filer-dg/agruppen/dg2/fengj/hbulbosum_genome/SNP_shortreads/ref_genome/FB19_011_3_pseudomolecules_v6_hap1.fasta'
bedtools makewindows -w 5000000 -g $ref.fai \
 | awk '{printf $0"\t"$1":"$2+1"-"$3; printf "\t%04d\n", NR}' > ${ref:t:r}_5Mb.bed

#run the SNP calling
interval='/filer-dg/agruppen/dg2/fengj/hbulbosum_genome/hb_gbs_analysis/GBSdata/ldhat/pollenseq/FB19_011_3_pseudomolecules_v6_hap1_5Mb.bed'
list='210504_bam_list.txt'
ref='/filer-dg/agruppen/dg2/fengj/hbulbosum_genome/SNP_shortreads/ref_genome/FB19_011_3_pseudomolecules_v6_hap1.fasta'
out='bcftools_210505'

mkdir $out

parallel -j 5 -a $interval --colsep '\t' ./run_bcftools_210504.zsh  $ref $list $out/$out '{4}' '{5}'

#convert VCF to long-format data.table
find bcftools_210505 | grep 'gz$' | sort | xargs zcat | ./normalize.awk | awk '$8 > 0' | awk 'length($3) == 1 && length($4) == 1' | bgzip > 210504_Hb_pollenseq_hap1_samtools_normalized.txt.gz
