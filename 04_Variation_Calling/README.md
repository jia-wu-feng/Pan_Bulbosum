# Pan_Bulbosum

Code for Article "A haplotype-resolved pangenome of the barley wild relative Hordeum bulbosum"

## Simulate long reads from haplotype-resolved genomes and call variants in each haplotype

````shell

for i in `cat genome_list.txt`;
do

echo $i
name=`basename $i`
echo $name
mkdir $name
cd $name

echo "#!/bin/bash

date


/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/final_genome/sim_long_reads/Badread/badread-runner.py simulate --reference $name.contig.fa --quantity 30x --error_model random \
         --qscore_model ideal --glitches 0,0,0 --junk_reads 0 --random_reads 0 --chimeras 0 \
         --identity 30,3 --length 40000,20000 --start_adapter_seq \"\" --end_adapter_seq \"\" \
  | gzip > $name.fastq.gz


pbmm2 align --preset HIFI --sort -j 26 -J 4 \
    --log-level INFO --log-file pbmm2.log --sample $name \
    /filer-dg/agruppen/dg7/fengj/hbulbosum_genome/final_genome/sim_long_reads/reference/FB19_011_3_pseudomolecules_v6_hap1.mmi $name.fastq.gz $name.bam 

samtools index -c $name.bam


mkdir tmp2
export TMPDIR=\"/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/final_genome/sim_long_reads/$name/tmp2\"
mkdir tmp
singularity exec -B /filer-dg/agruppen --pwd $PWD \
  /filer-dg/agruppen/dg7/fengj/hbulbosum_genome/final_genome/sim_long_reads/deepvariant_1.6.0.sif \
  run_deepvariant \
  --model_type=PACBIO \
  --ref=/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/final_genome/sim_long_reads/reference/FB19_011_3_pseudomolecules_v6_hap1.fasta \
  --reads=/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/final_genome/sim_long_reads/$name/$name.bam \
  --output_vcf=/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/final_genome/sim_long_reads/$name/$name.vcf.gz \
  --output_gvcf=/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/final_genome/sim_long_reads/$name/$name.g.vcf.gz \
  --intermediate_results_dir=/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/final_genome/sim_long_reads/$name/tmp \
  --make_examples_extra_args=vsc_min_count_indels=15 \
  --num_shards=40

date " > mapping.zsh

sbatch mapping.zsh

cd /filer-dg/agruppen/dg7/fengj/hbulbosum_genome/final_genome/sim_long_reads

done
````

## Merging and filtering variants

````shell

singularity exec -B /filer-dg/agruppen --pwd $PWD \
 /filer-dg/agruppen/dg7/fengj/hbulbosum_genome/final_genome/sim_long_reads/glnexus_v1.4.1.sif \
 /usr/local/bin/glnexus_cli \
 --config DeepVariant \
 --list /filer-dg/agruppen/dg7/fengj/hbulbosum_genome/final_genome/sim_long_reads/vcf_conbinant/gvcf_list.txt \
 --threads 8 -m 650 >  hbulbosum.deepv.filter.bcf

bcftools view hbulbosum.deepv.filter.bcf | bgzip -@ 8 -c  > hbulbosum.deepv.filter.vcf.gz

vcf='vcf_list.txt'
out='HB_SNP'
./recall_vcf.zsh --vcf $vcf --out $out --minmaf 0 --mindp 30 --minhomn 1 --dphom 3 --dphet 20 --minhomp 0 --threads 10 --tol 0.5 --minpresent 0 --minqual 20 > $out.out 2> $out.err

vcf='vcf_list.txt'
out='HB_Indel'
./recall_vcf_indel.zsh --vcf $vcf --out $out --minmaf 0 --mindp 30 --minhomn 1 --dphom 3 --dphet 20 --minhomp 0 --threads 10 --tol 0.5 --minpresent 0 --minqual 20 > $out.out 2> $out.err

````

## Call structural variations

````shell

for i in `cat genome_list.txt`;
do

echo $i
name=`basename $i`
nameb=`basename $i .fasta`
mkdir $name
echo $name

cd $name

ln -s /filer-dg/agruppen/dg7/fengj/hbulbosum_genome/final_genome/sim_long_reads/$name/$name.bam .

ln -s /filer-dg/agruppen/dg7/fengj/hbulbosum_genome/final_genome/sim_long_reads/$name/$name.bam.csi .


echo "#!/bin/bash

svim alignment --min_sv_size 50 $nameb $name.bam /filer-dg/agruppen/dg6/fengj/hbulbosum_genome/SV/genome/FB19_011_3_pseudomolecules_v6_hap1.fasta --insertion_sequences

" > svim.zsh

sbatch svim.zsh

cd /filer-dg/agruppen/dg6/fengj/hbulbosum_genome/SV 


done
````

## Merging and filtering structural variations


````shell

for i in `ls *.variants.vcf`; do echo $i ; name=`basename $i`; /filer-dg/agruppen/dg7/fengj/software/SURVIVOR-master/Debug/SURVIVOR filter $i NA 50 20000 0 5 ${i}.filter.vcf & ; done

ls *.filter.vcf > file.list

/filer-dg/agruppen/dg7/fengj/software/SURVIVOR-master/Debug/SURVIVOR merge file.list 1000 1 1 0 0 50 HB.SV.merge.vcf 

````
