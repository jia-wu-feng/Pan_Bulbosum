#!/bin/bash

############################### the parameters which users can reset ##########################################

# the absolute path of index file (.fai file) of the reference genome
fasta_index_file="/filer-dg/agruppen/dg2/fengj/hbulbosum_genome/IntroBlocker/FB19_011_3_pseudomolecules_v6_hap1.fasta.fai"

# chromosomes to analyze
chr_list_file="/filer-dg/agruppen/dg2/fengj/hbulbosum_genome/IntroBlocker/all_chr_ac0.65_5M/chr_list.txt"

# the absolute path of the working directory
working_dir="/filer-dg/agruppen/dg2/fengj/hbulbosum_genome/IntroBlocker/all_chr_ac0.65_5M"

# the absolute path of the IntroBlocker software directory
script_dir="/filer-dg/agruppen/dg2/fengj/hbulbosum_genome/IntroBlocker/all_chr_ac0.65_5M"

# the AHG bin size (in Mb)
bin_size=5

# the genomic variants density threshold used in the initial grouping (# of variants per Mb)
cluster_thr=8500

# the genomic variants density threshold used in the Bayesian smoothing (# of variants per Mb)
smooth_thr=1000

# the taxa_bam_file file contains taxa names and absolte paths of corresponding bam files of each chromosome
taxa_bam_file="/filer-dg/agruppen/dg2/fengj/hbulbosum_genome/IntroBlocker/all_chr_ac0.65_5M/taxa_bam.txt"

# the chr_vcf_file file contains absolte paths of vcf files of each chromosome
chr_vcf_file="/filer-dg/agruppen/dg2/fengj/hbulbosum_genome/IntroBlocker/all_chr_ac0.65_5M/chr_vcf.txt"

# the taxa_order_file file specifies the group of each taxa
taxa_order_file="/filer-dg/agruppen/dg2/fengj/hbulbosum_genome/IntroBlocker/all_chr_ac0.65_5M/taxa_order.txt"

# the mode of IntroBlocker (could be "un-supervised", "semi-supervised" and "supervised")
mode="semi-supervised"

############################## end of resetting parameters ##################################################
