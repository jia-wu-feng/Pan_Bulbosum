#!/bin/zsh

projectdir="/filer-dg/agruppen/dg2/fengj/hbulbosum_genome/genome_assembly/A40"

#In silico digest
#digest script is here: https://bitbucket.org/tritexassembly/tritexassembly.bitbucket.io/src/master/shell/
#Needs EMBOSS restrict, bedtools and rebase (http://rebase.neb.com), set paths in header

digest='/filer-dg/agruppen/seq_shared/mascher/code_repositories/tritexassembly.bitbucket.io/shell/digest_emboss.zsh'
ref='/filer-dg/agruppen/dg2/fengj/hbulbosum_genome/genome_assembly/A40/A40_hifiasm.fasta'
$digest --ref $ref --enzyme 'MboI' --sitelen 4 --minlen 30

##mapping script is here:
##https://bitbucket.org/tritexassembly/tritexassembly.bitbucket.io/src/master/shell/
##needs several standard tools, see header
##a free single-threaded version of Novosort is available from here:
##http://www.novocraft.com

map='/filer-dg/agruppen/seq_shared/mascher/code_repositories/tritexassembly.bitbucket.io/shell/run_hic_mapping.zsh'
bed=${ref:r}_MboI_fragments_30bp.bed

$map --threads 50 --mem '200G' --linker "GATCGATC" --ref $ref --bed $bed --tmp $TMPDIR /filer-dg/agruppen/dg2/fengj/hbulbosum_genome/genome_assembly/A40/hic/hicreads

##Use AAGCTAGCTT for HindIII


