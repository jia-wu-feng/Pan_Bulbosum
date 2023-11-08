#!/bin/zsh

projectdir="/filer-dg/agruppen/dg2/fengj/hbulbosum_genome/genome_assembly/A40"

#Index construction

fasta='/filer-dg/agruppen/dg2/fengj/hbulbosum_genome/genome_assembly/A40/A40_hifiasm.fasta'
gmap_build='/opt/Bio/gmap/2019-09-12/bin/gmap_build'
dir=$fasta:h'/gmap_index'
name=${fasta:t:r}
$gmap_build $fasta -D $dir -d $name > ${name}_gmap_build.out  2> ${name}_gmap_build.err

#Alignment
gmapl='/opt/Bio/gmap/2019-09-12/bin/gmapl'
gmap_build='/opt/Bio/gmap/2019-09-12/bin/gmap_build'
dir=$fasta:h'/gmap_index'
name=${fasta:t:r}
threads=30
query='/filer-dg/agruppen/seq_shared/mascher/morex_v3_annotation_200818/Hv_Morex.pgsb.Jul2020.HC.cds_longest.fa'
prefix="${name}_MorexV3_HC"

$gmapl -d $name -D $dir -t $threads -f 2 $query > ${prefix}.gff 2> ${prefix}.err  
grep mRNA $prefix.gff | cut -f 1,4,5,9 | tr ';' '\t' | cut -f 1-3,5,7-8 \
   | sed -e 's/coverage=//' -e 's/identity=//' -e 's/Name=//' > ${prefix}_table.txt 



