# Pan_Bulbosum

Code for Article "A haplotype-resolved pangenome of the barley wild relative Hordeum bulbosum"

## Diploid genome assembly (example accession: FB20-005-1)

### Step 1 

#### Genome size evaluation 

````shell
zcat ../hifi_reads/\*.fastq.gz | jellyfish count /dev/fd/0 -C -o FB20_005_1_k71 -m 71 -t 20 -s 5G

jellyfish histo -h 3000000 -o FB_19_028_3_k71.histo FB20_005_1_k71
````

````R
.libPaths(c("/home/mascher/tmp/3.5.1", "/filer-dg/agruppen/seq_shared/mascher/Rlib/3.5.1", "/opt/Bio/R_LIBS/3.5.1"))

library("findGSE")

findGSE(histo="FB20_005_1_k71.histo", sizek=71, outdir="FB20_005_1_he_k71", exp_hom=80)
````

#### Assembly using HiFiasm

Hifiasm.zsh

Convert gfa to fasta and calculate N50

````shell
gfatools gfa2fa FB20_005_1_hifiasm.p_utg.gfa | seqtk rename - 'contig_' > FB20_005_1_hifiasm.p_utg.fasta

samtools faidx FB20_005_1_hifiasm.p_utg.fasta

/filer-dg/agruppen/seq_shared/mascher/code_repositories/tritexassembly.bitbucket.io/shell/n50 FB20_005_1_hifiasm.p_utg.fasta.fai > FB20_005_1_hifiasm.p_utg.fasta.n50 


````

#### Align Hi-C data to assembled unitigs

Hi-C_map.zsh

#### Use the barley high-confidence gene as the maker to build a guide map

Guide_map.zsh

#### create assembly object 

Create_assembly_object.R


### Step 2

Phasing.R

### Step 3

Output_assembly.R

Output unassembled sequence

````R
source('/filer-dg/agruppen/seq_shared/mascher/code_repositories/triticeae.bitbucket.io/R/pseudomolecule_construction.R')

readRDS('FB20_005_1_hic_map_v2_hap1.Rds') -> hic_map_v2_hap1
readRDS('FB20_005_1_hic_map_v2_hap2.Rds') -> hic_map_v2_hap2



hic_map_v2_hap1$agp[!is.na(chr)]$scaffold->h1
hic_map_v2_hap2$agp[!is.na(chr)]$scaffold->h2


hic_map_v2_hap1$agp[,.(scaffold,scaffold_length)][scaffold!="gap"] -> chrh

chrh[,hap:=0]
chrh[scaffold %in% h1, hap:=hap+10]
chrh[scaffold %in% h2, hap:=hap+2]


write.table(chrh,"FB20_005_1_contig_all.txt",quote=F,row.names=F,sep="\t")

sum(hic_map_v2_hap1$agp[!is.na(chr)]$scaffold_length)
sum(hic_map_v2_hap2$agp[!is.na(chr)]$scaffold_length)


sum(chrh[hap==0]$scaffold_length)

````


````shell
less FB20_005_1_contig_all.txt | awk '{if($3==0) print $1}' > uncontig

seqkit grep -f uncontig FB20_005_1_pseudomolecules_v2_hap1/*_assembly_v2.fasta > FB20_005_1_unanchor.fasta


````
