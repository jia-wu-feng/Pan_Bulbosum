#https://bitbucket.org/tritexassembly/tritexassembly.bitbucket.io/src/master/R/
source('/filer-dg/agruppen/seq_shared/mascher/code_repositories/triticeae.bitbucket.io/R/pseudomolecule_construction.R')

#read chromosome lengths of Morex V3
fread('/filer-dg/agruppen/seq_shared/mascher/morexV3_pseudomolecules_200421/200416_MorexV3_pseudomolecules.fasta.fai', sel=1:2, col.names=c("chr", "len"))->fai

#read centromere positions
fread('/filer-dg/agruppen/seq_shared/mascher/morexV3_pseudomolecules_200421/MorexV3_centromere_positions.tsv')->cen
setnames(cen, c("chr", "cen_pos"))

#read Hi-C mapping output (list of read pairs connecting RE fragments)
fread('/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/FB20_005_1/hic/hicreads/hicreads_fragment_pairs.tsv.gz')->f
setnames(f, c("ctg1", "pos1", "ctg2", "pos2"))

#read contig lengths of assebmbly
faif <- '/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/FB20_005_1/FB20_005_1_hifiasm.p_utg.fasta.fai'
fread(faif, sel=1:2, col.names=c("contig", "contig_length"))->fb

#read GMAP alignments of Morex V3 genes
fread('./FB20_005_1_hifiasm.p_utg_MorexV3_HC_table.txt')->x
setnames(x, c("contig", "start", "end", "transcript", "alnlen", "id"))

#read position of genes in Morex assembly
fread('/filer-dg/agruppen/seq_shared/mascher/morex_v3_annotation_200818/Hv_Morex.pgsb.Jul2020.HC_mrna_pos.txt')->p
setnames(p, c("chr", "chr_start", "chr_end", "transcript"))

#merge mapping with gene positions
p[x, on="transcript"] -> px
#filter for high-quality alignments
px[alnlen >= 80 & id >= 90]->a
#only consider genes with up to 2 alignments (one to each haplotype)
px[a[, .N, key=transcript][N <= 2]$transcript, on="transcript"]->aa

#keep contigs with at least 4 aligned genes, at least 75 % of which come from the same chromosome
aa[, .N, key=.(chr, contig)][, p := N/sum(N), by=contig][order(-p)][!duplicated(contig)] -> cc
cc[N >= 4 & p >= 0.75] -> cc

#merge with contig lengths
fb[cc, on="contig"] -> cc

#check proportion of anchored sequence
sum(cc$contig_length)/1e6

#calculate approximate chromosomal locations
aa[cc[, .(contig, chr)], on=c("chr", "contig")][, .(pos=median(as.numeric(chr_start)), pos_mad=mad(as.numeric(chr_start))), key=contig][cc, on="contig"] -> cc
cc[, mr := pos_mad / contig_length]
setorder(cc, chr, pos) 
cc[, chr_idx := 1:.N, by=chr]
cc[, agp_pos := c(0, cumsum(contig_length[-.N])), by=chr]

#read coverage information and convert to correct assembly coordinates
fread(cmd="grep '^S' ../FB20_005_1_hifiasm.p_utg.noseq.gfa | cut -f 2,5 | tr ':' '\\t' | cut -f 1,4", head=F) -> cov
setnames(cov, c("contig", "cc"))
cov[, contig := sub('l$', '', sub("utg0*", "contig_", contig))]
cc[cov, on="contig"] ->  cc

#### Get the distribution of coverages, infer 1X and 2X boundaries

pdf("coverage_chr.pdf", height=10, width=10)
pcc1<-cc[chr!="<NA>"]
plot(density(pcc1$cc),xlim=c(0,100))
dev.off()

##
# Import data and create assembly object
##

#read MorexV3 gene-based guide map ("pseudo-POPSEQ")
readRDS('/filer-dg/agruppen/seq_shared/mascher/hordeum_bulbosum_fb19_011_ccs_assembly_201111/Hb_FB19_011_hifiasm_201217/210129_MorexV3_HCgenes_pseudopopseq.Rds') -> pg

#read unitig lengths
f <- '/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/FB20_005_1/FB20_005_1_hifiasm.p_utg.fasta.fai'
fread(f, head=F, select=1:2, col.names=c("scaffold", "length"))->fai

#read MorexV3 gene alignment (GMAP output) and merge with pseudo-POPSEQ table
fread('/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/FB20_005_1/hic/FB20_005_1_hifiasm.p_utg_MorexV3_HC_table.txt')->x
setnames(x, c("contig", "start", "end", "transcript", "alnlen", "id"))
fai[, .(scaffold, scaffold_length=length)][x[, .(css_contig=transcript, scaffold=contig, pos=start)], on="scaffold"] -> aln
pg[aln, on="css_contig"] -> aln

#read Hi-C map
dir <- '/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/FB20_005_1/hic/hicreads'
fread(paste('find', dir, '| grep "fragment_pairs.tsv.gz$" | xargs zcat'),
      header=F, col.names=c("scaffold1", "pos1", "scaffold2", "pos2"))->fpairs

#initialize assembly object, anchor to guide map and calculate physical coverage
init_assembly(fai=fai, cssaln=aln, fpairs=fpairs) -> assembly
anchor_scaffolds(assembly = assembly, popseq=pg, species="barley") -> assembly
add_hic_cov(assembly, binsize=1e4, binsize2=1e6, minNbin=50, innerDist=3e5, cores=20)->assembly

#save uncorrected assembly object
saveRDS(assembly, file="FB20_005_1_assembly.Rds")

##
# Check for and correct chimeras 
##

readRDS('FB20_005_1_assembly.Rds') -> assembly

#make diagnostic plot for one contig


#plot all > 1e6 contig & mri <= -2.5 
assembly$info[length >= 1e6& mri <= -2.5, .(scaffold, length)][order(-length)] -> ss

i=1
ss[i]$scaffold -> s 
assembly$cov[s, on='scaffold'] -> b

for (ni in 2:length(ss$scaffold)) {
i=ni
ss[i]$scaffold -> s 
rbind(b, assembly$cov[s, on='scaffold'])->b
} 

plot_chimeras(assembly=assembly, scaffolds=b,  species="barley", refname="Morex HC genes", autobreaks=F, mbscale=1,file="FB20_005_1_assembly_all_contig.pdf", cores=50)

#set break point coordinates

i=5
ss[i]$scaffold -> s 
assembly$cov[s, on='scaffold'][bin >= 10e6 & bin <= 12e6][order(r)][1, .(scaffold, bin)] -> b

i=44
ss[i]$scaffold -> s 
rbind(b, assembly$cov[s, on='scaffold'][bin >= 3e6 & bin <= 4e6][order(r)][1, .(scaffold, bin)])->b


setnames(b, "bin", "br")
plot_chimeras(assembly=assembly, scaffolds=b, br=b, species="barley", refname="MorexV3 genes",  mbscale=1,
         file="FB20_005_1_assembly_chimeras_final.pdf", cores=30)

#implement the correction
break_scaffolds(b, assembly, prefix="contig_corrected_v1_", slop=1e4, cores=30, species="barley") -> assembly_v2

#save the object
saveRDS(assembly_v2, file="FB20_005_1_assembly_v2.Rds")
