#https://bitbucket.org/tritexassembly/tritexassembly.bitbucket.io/src/master/R/
source('/filer-dg/agruppen/seq_shared/mascher/code_repositories/triticeae.bitbucket.io/R/pseudomolecule_construction.R')

#read chromosome lengths of Morex V3
fread('/filer-dg/agruppen/seq_shared/mascher/morexV3_pseudomolecules_200421/200416_MorexV3_pseudomolecules.fasta.fai', sel=1:2, col.names=c("chr", "len"))->fai

#read centromere positions
fread('/filer-dg/agruppen/seq_shared/mascher/morexV3_pseudomolecules_200421/MorexV3_centromere_positions.tsv')->cen
setnames(cen, c("chr", "cen_pos"))

readRDS('A40_assembly_v2.Rds') -> assembly_v2


chrNames <- function(agp=F, species="wheat") {
 if(species == "wheat"){
  data.table(alphachr=apply(expand.grid(1:7, c("A", "B", "D"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:21)->z
 } else if (species == "barley"){
  data.table(alphachr=apply(expand.grid(1:7, "H", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "rye"){
  data.table(alphachr=apply(expand.grid(1:7, "R", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "lolium"){
  data.table(alphachr=as.character(1:7), chr=1:7)->z
 }
 else if (species == "maize"){
  data.table(alphachr=as.character(1:10), chr=1:10)->z
 }
 else if (species == "sharonensis"){
  data.table(alphachr=apply(expand.grid(1:7, "S", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "oats"){
  data.table(alphachr=sub(" ", "", apply(expand.grid(1:21, "M", stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:21)->z
 }
 else if (species == "oats_new"){
  data.table(alphachr=apply(expand.grid(1:7, c("A", "C", "D"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:21)->z
 }
 else if (species == "avena_barbata"){
  data.table(alphachr=apply(expand.grid(1:7, c("A", "B"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:14)->z
 }
 else if (species == "hordeum_bulbosum"){
  data.table(alphachr=sort(apply(expand.grid(1:7, c("H_"), c(1:2), stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:14)->z
 }
 else if (species == "hordeum_bulbosum_4x"){
  data.table(alphachr=sort(apply(expand.grid(1:7, c("H_"), c(1:4), stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:28)->z
 }
 if(agp){
  rbind(z, data.table(alphachr="Un", chr=0))[, agp_chr := paste0("chr", alphachr)]->z
 }
 z[]
}

wheatchr <- function(agp=F, species="wheat") {
 if(species == "wheat"){
  data.table(alphachr=apply(expand.grid(1:7, c("A", "B", "D"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:21)->z
 } else if (species == "barley"){
  data.table(alphachr=apply(expand.grid(1:7, "H", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "rye"){
  data.table(alphachr=apply(expand.grid(1:7, "R", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "lolium"){
  data.table(alphachr=as.character(1:7), chr=1:7)->z
 }
 else if (species == "maize"){
  data.table(alphachr=as.character(1:10), chr=1:10)->z
 }
 else if (species == "sharonensis"){
  data.table(alphachr=apply(expand.grid(1:7, "S", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "oats"){
  data.table(alphachr=sub(" ", "", apply(expand.grid(1:21, "M", stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:21)->z
 }
 else if (species == "oats_new"){
  data.table(alphachr=apply(expand.grid(1:7, c("A", "C", "D"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:21)->z
 }
 else if (species == "avena_barbata"){
  data.table(alphachr=apply(expand.grid(1:7, c("A", "B"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:14)->z
 }
 else if (species == "hordeum_bulbosum"){
  data.table(alphachr=sort(apply(expand.grid(1:7, c("H_"), c(1:2), stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:14)->z
 }
 else if (species == "hordeum_bulbosum_4x"){
  data.table(alphachr=sort(apply(expand.grid(1:7, c("H_"), c(1:4), stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:28)->z
 }
 if(agp){
  rbind(z, data.table(alphachr="Un", chr=0))[, agp_chr := paste0("chr", alphachr)]->z
 }
 z[]
}



#Function to combine the Hi-C maps from two haplotypes
combine_hic <- function(hh,hap1, hap2, hap3, hap4, assembly, species="hordeum_bulbosum_4x"){
 assembly_v2 <- assembly
 
 hic_map_v1_hap1 <- hap1
 hic_map_v1_hap2 <- hap2
 hic_map_v1_hap3 <- hap3
 hic_map_v1_hap4 <- hap4

 hic_map_v1_hap1$agp[agp_chr != "chrUn"] -> a1
 hic_map_v1_hap2$agp[agp_chr != "chrUn"] -> a2
 hic_map_v1_hap3$agp[agp_chr != "chrUn"] -> a3
 hic_map_v1_hap4$agp[agp_chr != "chrUn"] -> a4
 
 a1[, agp_chr := paste0(agp_chr, "_1")]
 a2[, agp_chr := paste0(agp_chr, "_2")]
 a3[, agp_chr := paste0(agp_chr, "_3")]
 a4[, agp_chr := paste0(agp_chr, "_4")]
 
 a1[, chr := NULL]
 a2[, chr := NULL]
 a3[, chr := NULL]
 a4[, chr := NULL]

 chrNames(agp=T, species)[, .(chr, agp_chr)][a1, on="agp_chr"] -> a1 ###### there for ask  
 chrNames(agp=T, species)[, .(chr, agp_chr)][a2, on="agp_chr"] -> a2
 chrNames(agp=T, species)[, .(chr, agp_chr)][a3, on="agp_chr"] -> a3
 chrNames(agp=T, species)[, .(chr, agp_chr)][a4, on="agp_chr"] -> a4

 c(a1[scaffold != "gap"]$scaffold, a2[scaffold != "gap"]$scaffold, a3[scaffold != "gap"]$scaffold, a4[scaffold != "gap"]$scaffold) -> s
 s[duplicated(s)] -> s


 a1[hh[grepl("1", hap) & hap>=10]$scaffold, on="scaffold", scaffold := paste0(scaffold, "_hap1")]
 a2[hh[grepl("2", hap) & hap>=10]$scaffold, on="scaffold", scaffold := paste0(scaffold, "_hap2")]
 a3[hh[grepl("3", hap) & hap>=10]$scaffold, on="scaffold", scaffold := paste0(scaffold, "_hap3")]
 a4[hh[grepl("4", hap) & hap>=10]$scaffold, on="scaffold", scaffold := paste0(scaffold, "_hap4")]
 rbind(a1, a2, a3, a4) -> a
 
 hic_map_v1_hap1$chrlen[!is.na(chr)][, .(agp_chr=paste0(agp_chr, "_1"), length, truechr)] -> l1
 hic_map_v1_hap2$chrlen[!is.na(chr)][, .(agp_chr=paste0(agp_chr, "_2"), length, truechr)] -> l2
 hic_map_v1_hap3$chrlen[!is.na(chr)][, .(agp_chr=paste0(agp_chr, "_3"), length, truechr)] -> l3
 hic_map_v1_hap4$chrlen[!is.na(chr)][, .(agp_chr=paste0(agp_chr, "_4"), length, truechr)] -> l4

 rbind(l1, l2, l3, l4) -> l
 l[, offset := cumsum(c(0, length[1:(.N-1)]))]
 l[, plot_offset := cumsum(c(0, length[1:(.N-1)]+1e8))]
 chrNames(agp=T, species)[l, on="agp_chr"] -> l

 copy(assembly_v2$info) -> ai
 ai[!s, on="scaffold"] -> u
 ai[hh[grepl("1", hap) & hap >= 10]$scaffold, on="scaffold"][, scaffold := paste0(scaffold, "_hap1")] -> i1
 ai[hh[grepl("2", hap) & hap >= 10]$scaffold, on="scaffold"][, scaffold := paste0(scaffold, "_hap2")] -> i2
 ai[hh[grepl("3", hap) & hap >= 10]$scaffold, on="scaffold"][, scaffold := paste0(scaffold, "_hap3")] -> i3
 ai[hh[grepl("4", hap) & hap >= 10]$scaffold, on="scaffold"][, scaffold := paste0(scaffold, "_hap4")] -> i4

 rbind(u, i1, i2, i3, i4) -> i

 assembly_v2$fpairs[, .(scaffold1, scaffold2, pos1, pos2)] -> f
 f[scaffold1 %in% hh[hap==1234]$scaffold, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.5, ifelse(runif(.N) > 0.5, "_hap1","_hap2"), ifelse(runif(.N) > 0.5, "_hap3","_hap4" )))] 
 f[scaffold1 %in% hh[hap==123]$scaffold, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.333, ifelse(runif(.N) > 0.5, "_hap1","_hap2"), "_hap3"))]
 f[scaffold1 %in% hh[hap==134]$scaffold, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.333, ifelse(runif(.N) > 0.5, "_hap1","_hap3"), "_hap4"))]
 f[scaffold1 %in% hh[hap==124]$scaffold, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.333, ifelse(runif(.N) > 0.5, "_hap1","_hap2"), "_hap4"))]
 f[scaffold1 %in% hh[hap==234]$scaffold, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.333, ifelse(runif(.N) > 0.5, "_hap2","_hap3"), "_hap4"))]
 f[scaffold1 %in% hh[hap==12]$scaffold, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.5, "_hap1", "_hap2"))]
 f[scaffold1 %in% hh[hap==13]$scaffold, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.5, "_hap1", "_hap3"))]
 f[scaffold1 %in% hh[hap==14]$scaffold, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.5, "_hap1", "_hap4"))]
 f[scaffold1 %in% hh[hap==23]$scaffold, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.5, "_hap2", "_hap3"))]
 f[scaffold1 %in% hh[hap==24]$scaffold, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.5, "_hap2", "_hap4"))]
 f[scaffold1 %in% hh[hap==34]$scaffold, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.5, "_hap3", "_hap4"))]
 

 f[scaffold2 %in% hh[hap==1234]$scaffold, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.5, ifelse(runif(.N) > 0.5, "_hap1","_hap2"), ifelse(runif(.N) > 0.5, "_hap3","_hap4" )))] 
 f[scaffold2 %in% hh[hap==123]$scaffold, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.333, ifelse(runif(.N) > 0.5, "_hap1","_hap2"), "_hap3"))]
 f[scaffold2 %in% hh[hap==134]$scaffold, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.333, ifelse(runif(.N) > 0.5, "_hap1","_hap3"), "_hap4"))]
 f[scaffold2 %in% hh[hap==124]$scaffold, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.333, ifelse(runif(.N) > 0.5, "_hap1","_hap2"), "_hap4"))]
 f[scaffold2 %in% hh[hap==234]$scaffold, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.333, ifelse(runif(.N) > 0.5, "_hap2","_hap3"), "_hap4"))]
 f[scaffold2 %in% hh[hap==12]$scaffold, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.5, "_hap1", "_hap2"))]
 f[scaffold2 %in% hh[hap==13]$scaffold, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.5, "_hap1", "_hap3"))]
 f[scaffold2 %in% hh[hap==14]$scaffold, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.5, "_hap1", "_hap4"))]
 f[scaffold2 %in% hh[hap==23]$scaffold, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.5, "_hap2", "_hap3"))]
 f[scaffold2 %in% hh[hap==24]$scaffold, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.5, "_hap2", "_hap4"))]
 f[scaffold2 %in% hh[hap==34]$scaffold, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.5, "_hap3", "_hap4"))]

 list(info=i, fpairs=f) -> assembly_hap
 list(agp=a, chrlen=l) -> hic_map
 list(assembly_hap=assembly_hap, hic_map=hic_map)
}


##
# Construct first Hi-C map
##

#read Hi-C fragment data
f <- '../A40_hifiasm_MboI_fragments_30bp.bed'
read_fragdata(info=assembly_v2$info, file=f)->frag_data

#read haplotype information
readRDS(file="A40_assembly_v2_Hv_guide+HiClift+hap1.Rds") -> hh

###for 4 hap 
##c(1,2,3,4,12,13,14,23,24,34,123,124,134,234,1234)

#Hi-C map for haplotype 1
hh[hap %in% c(1,12,13,14,123,124,134,1234), .(scaffold, length, chr=as.integer(substr(chr, 4, 4)), cM = pos/1e6)] -> hic_info
frag_data$info[, .(scaffold, nfrag)][hic_info, on="scaffold"] -> hic_info

hic_map(info=hic_info, assembly=assembly_v2, frags=frag_data$bed, species="barley", ncores=21,
    min_nfrag_scaffold=30, max_cM_dist = 450,
    binsize=2e5, min_nfrag_bin=10, gap_size=100)->hic_map_v1_hap1

saveRDS(hic_map_v1_hap1, file="A40_hifiasm_map_v1_hap1.Rds")

#Hi-C map for haplotype 2
hh[hap %in% c(2,12,23,24,123,124,234,1234), .(scaffold, length, chr=as.integer(substr(chr, 4, 4)), cM = pos/1e6)] -> hic_info
frag_data$info[, .(scaffold, nfrag)][hic_info, on="scaffold"] -> hic_info

hic_map(info=hic_info, assembly=assembly_v2, frags=frag_data$bed, species="barley", ncores=21,
    min_nfrag_scaffold=30, max_cM_dist = 450,
    binsize=2e5, min_nfrag_bin=10, gap_size=100)->hic_map_v1_hap2

saveRDS(hic_map_v1_hap2, file="A40_hifiasm_map_v1_hap2.Rds")


#Hi-C map for haplotype 3
hh[hap %in% c(3,13,23,34,123,134,234,1234), .(scaffold, length, chr=as.integer(substr(chr, 4, 4)), cM = pos/1e6)] -> hic_info
frag_data$info[, .(scaffold, nfrag)][hic_info, on="scaffold"] -> hic_info

hic_map(info=hic_info, assembly=assembly_v2, frags=frag_data$bed, species="barley", ncores=21,
    min_nfrag_scaffold=30, max_cM_dist = 450,
    binsize=2e5, min_nfrag_bin=10, gap_size=100)->hic_map_v1_hap3

saveRDS(hic_map_v1_hap3, file="A40_hifiasm_map_v1_hap3.Rds")


#Hi-C map for haplotype 4
hh[hap %in% c(4,14,24,34,124,134,234,1234), .(scaffold, length, chr=as.integer(substr(chr, 4, 4)), cM = pos/1e6)] -> hic_info
frag_data$info[, .(scaffold, nfrag)][hic_info, on="scaffold"] -> hic_info

hic_map(info=hic_info, assembly=assembly_v2, frags=frag_data$bed, species="barley", ncores=21,
    min_nfrag_scaffold=30, max_cM_dist = 450,
    binsize=2e5, min_nfrag_bin=10, gap_size=100)->hic_map_v1_hap4

saveRDS(hic_map_v1_hap4, file="A40_hifiasm_map_v1_hap4.Rds")




#create Hi-C plots for haplotype 1
snuc <- '../A40_hifiasm_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="A40_hifiasm_map_v1_hap1.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v1_hap1_l

#create Hi-C plots for haplotype 2
snuc <- '../A40_hifiasm_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="A40_hifiasm_map_v1_hap2.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v1_hap2_l

#create Hi-C plots for haplotype 3
snuc <- '../A40_hifiasm_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="A40_hifiasm_map_v1_hap3.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v1_hap3_l

#create Hi-C plots for haplotype 4
snuc <- '../A40_hifiasm_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="A40_hifiasm_map_v1_hap4.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v1_hap4_l


#combine both Hi-C maps
readRDS("A40_hifiasm_map_v1_hap1.Rds")->hic_map_v1_hap1
readRDS("A40_hifiasm_map_v1_hap2.Rds")->hic_map_v1_hap2
readRDS("A40_hifiasm_map_v1_hap3.Rds")->hic_map_v1_hap3
readRDS("A40_hifiasm_map_v1_hap4.Rds")->hic_map_v1_hap4


combine_hic(hh,hic_map_v1_hap1, hic_map_v1_hap2, hic_map_v1_hap3, hic_map_v1_hap4,assembly_v2) -> ch
ch$hic_map -> hic_map_v1
ch$assembly_hap -> assembly_v2_hap

#create inter-chromosomal plot for combined hap1+hap2 Hi-C haplotypes
nuc <- '../A40_hifiasm_MboI_fragments_30bp_split.nuc.txt'
add_psmol_fpairs(assembly=assembly_v2_hap, hic_map=hic_map_v1, nucfile=nuc)->hic_map_v1_l

bin_hic_step(hic=hic_map_v1_l$links, frags=hic_map_v1_l$frags, binsize=1e6, chrlen=hic_map_v1_l$chrlen, cores=14)->hic_map_v1_l$hic_1Mb

f <- "A40_hifiasm_map_v1_interchromosomal.png"
interchromosomal_matrix_plot(hic_map=hic_map_v1_l, file=f, species="hordeum_bulbosum_4x")

#write Hi-C map for editing
write_hic_map(rds="A40_hifiasm_map_v1_hap1.Rds", file="A40_hifiasm_map_v1_hap1.xlsx", species="barley")
write_hic_map(rds="A40_hifiasm_map_v1_hap2.Rds", file="A40_hifiasm_map_v1_hap2.xlsx", species="barley")
write_hic_map(rds="A40_hifiasm_map_v1_hap3.Rds", file="A40_hifiasm_map_v1_hap3.xlsx", species="barley")
write_hic_map(rds="A40_hifiasm_map_v1_hap4.Rds", file="A40_hifiasm_map_v1_hap4.xlsx", species="barley")

###### write fasta

readRDS('A40_assembly_v2.Rds') -> assembly_v2
readRDS('A40_hifiasm_map_v1_hap1.Rds') -> hic_map_v1_hap1
readRDS('A40_hifiasm_map_v1_hap2.Rds') -> hic_map_v1_hap2
readRDS('A40_hifiasm_map_v1_hap3.Rds') -> hic_map_v1_hap3
readRDS('A40_hifiasm_map_v1_hap4.Rds') -> hic_map_v1_hap4


fasta <- '/filer-dg/agruppen/dg2/fengj/hbulbosum_genome/genome_assembly/A40/A40_hifiasm.fasta'
sink("A40_pseudomolecules_v1_hap1.log")
compile_psmol(fasta=fasta, output="A40_pseudomolecules_v1_hap1", hic_map=hic_map_v1_hap1, assembly=assembly_v2, cores=30)
sink()

sink("A40_pseudomolecules_v1_hap2.log")
compile_psmol(fasta=fasta, output="A40_pseudomolecules_v1_hap2", hic_map=hic_map_v1_hap2, assembly=assembly_v2, cores=30)
sink()

sink("A40_pseudomolecules_v1_hap3.log")
compile_psmol(fasta=fasta, output="A40_pseudomolecules_v1_hap3", hic_map=hic_map_v1_hap3, assembly=assembly_v2, cores=30)
sink()

sink("A40_pseudomolecules_v1_hap4.log")
compile_psmol(fasta=fasta, output="A40_pseudomolecules_v1_hap4", hic_map=hic_map_v1_hap4, assembly=assembly_v2, cores=30)
sink()



### Perform mannel check for each haplotype



#import Hi-C maps after editing, v1 -> v2
read_hic_map(rds="A40_hifiasm_map_v1_hap1.Rds", file="A40_hifiasm_map_v1_hap1_edit.xlsx") -> nmap
hic_map(species="barley", agp_only=T, map=nmap)->hic_map_v2_hap1
saveRDS(hic_map_v2_hap1, file="A40_hifiasm_map_v2_hap1.Rds")

read_hic_map(rds="A40_hifiasm_map_v1_hap2.Rds", file="A40_hifiasm_map_v1_hap2_edit.xlsx") -> nmap
hic_map(species="barley", agp_only=T, map=nmap)->hic_map_v2_hap2
saveRDS(hic_map_v2_hap2, file="A40_hifiasm_map_v2_hap2.Rds")

read_hic_map(rds="A40_hifiasm_map_v1_hap3.Rds", file="A40_hifiasm_map_v1_hap3_edit.xlsx") -> nmap
hic_map(species="barley", agp_only=T, map=nmap)->hic_map_v2_hap3
saveRDS(hic_map_v2_hap3, file="A40_hifiasm_map_v2_hap3.Rds")

read_hic_map(rds="A40_hifiasm_map_v1_hap4.Rds", file="A40_hifiasm_map_v1_hap4_edit.xlsx") -> nmap
hic_map(species="barley", agp_only=T, map=nmap)->hic_map_v2_hap4
saveRDS(hic_map_v2_hap4, file="A40_hifiasm_map_v2_hap4.Rds")


#proceed with Hi-C plots and repeat edit cycle if need be


#create Hi-C plots for haplotype 1
snuc <- '../A40_hifiasm_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="A40_hifiasm_map_v2_hap1.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v2_hap1_l

#create Hi-C plots for haplotype 2
snuc <- '../A40_hifiasm_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="A40_hifiasm_map_v2_hap2.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v1_hap2_l

#create Hi-C plots for haplotype 3
snuc <- '../A40_hifiasm_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="A40_hifiasm_map_v2_hap3.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v1_hap3_l

#create Hi-C plots for haplotype 4
snuc <- '../A40_hifiasm_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="A40_hifiasm_map_v2_hap4.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v1_hap4_l


#readRDS('A40_assembly_v2.Rds') -> assembly_v2
readRDS('A40_hifiasm_map_v2_hap1.Rds') -> hic_map_v2_hap1
readRDS('A40_hifiasm_map_v2_hap2.Rds') -> hic_map_v2_hap2
readRDS('A40_hifiasm_map_v2_hap3.Rds') -> hic_map_v2_hap3
readRDS('A40_hifiasm_map_v2_hap4.Rds') -> hic_map_v2_hap4


fasta <- '/filer-dg/agruppen/dg2/fengj/hbulbosum_genome/genome_assembly/A40/A40_hifiasm.fasta'

sink("A40_pseudomolecules_v2_hap1.log")
compile_psmol(fasta=fasta, output="A40_pseudomolecules_v2_hap1", hic_map=hic_map_v2_hap1, assembly=assembly_v2, cores=30)
sink()

sink("A40_pseudomolecules_v2_hap2.log")
compile_psmol(fasta=fasta, output="A40_pseudomolecules_v2_hap2", hic_map=hic_map_v2_hap2, assembly=assembly_v2, cores=30)
sink()

sink("A40_pseudomolecules_v2_hap3.log")
compile_psmol(fasta=fasta, output="A40_pseudomolecules_v2_hap3", hic_map=hic_map_v2_hap3, assembly=assembly_v2, cores=30)
sink()

sink("A40_pseudomolecules_v2_hap4.log")
compile_psmol(fasta=fasta, output="A40_pseudomolecules_v2_hap4", hic_map=hic_map_v2_hap4, assembly=assembly_v2, cores=30)
sink()


readRDS(file="A40_assembly_v2_Hv_guide+HiClift+hap1.Rds") -> hh

combine_hic(hh,hic_map_v2_hap1, hic_map_v2_hap2, hic_map_v2_hap3, hic_map_v2_hap4,assembly_v2) -> ch2
ch2$hic_map -> hic_map_v2
ch2$assembly_hap -> assembly_v2_hap

#create inter-chromosomal plot for combined hap1+hap2 Hi-C haplotypes
nuc <- '../A40_hifiasm_MboI_fragments_30bp_split.nuc.txt'
add_psmol_fpairs(assembly=assembly_v2_hap, hic_map=hic_map_v2, nucfile=nuc)->hic_map_v2_l

bin_hic_step(hic=hic_map_v2_l$links, frags=hic_map_v2_l$frags, binsize=1e6, chrlen=hic_map_v2_l$chrlen, cores=21)->hic_map_v2_l$hic_1Mb

f <- "A40_hifiasm_map_v2_interchromosomal.png"
interchromosomal_matrix_plot(hic_map=hic_map_v2_l, file=f, species="hordeum_bulbosum_4x")


