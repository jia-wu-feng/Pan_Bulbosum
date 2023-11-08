#https://bitbucket.org/tritexassembly/tritexassembly.bitbucket.io/src/master/R/
source('/filer-dg/agruppen/seq_shared/mascher/code_repositories/triticeae.bitbucket.io/R/pseudomolecule_construction.R')

#Function to combine the Hi-C maps from two haplotypes
combine_hic <- function(hap1, hap2, assembly, species="hordeum_bulbosum"){
 assembly_v2 <- assembly
 hic_map_v1_hap1 <- hap1
 hic_map_v1_hap2 <- hap2
 hic_map_v1_hap1$agp[agp_chr != "chrUn"] -> a1
 hic_map_v1_hap2$agp[agp_chr != "chrUn"] -> a2
 a1[, agp_chr := paste0(agp_chr, "_1")]
 a2[, agp_chr := paste0(agp_chr, "_2")]
 a1[, chr := NULL]
 a2[, chr := NULL]
 chrNames(agp=T, species)[, .(chr, agp_chr)][a1, on="agp_chr"] -> a1
 chrNames(agp=T, species)[, .(chr, agp_chr)][a2, on="agp_chr"] -> a2

 c(a1[scaffold != "gap"]$scaffold, a2[scaffold != "gap"]$scaffold) -> s
 s[duplicated(s)] -> s

 a1[s, on="scaffold", scaffold := paste0(scaffold, "_hap1")]
 a2[s, on="scaffold", scaffold := paste0(scaffold, "_hap2")]
 rbind(a1, a2) -> a
 
 hic_map_v1_hap1$chrlen[!is.na(chr)][, .(agp_chr=paste0(agp_chr, "_1"), length, truechr)] -> l1
 hic_map_v1_hap2$chrlen[!is.na(chr)][, .(agp_chr=paste0(agp_chr, "_2"), length, truechr)] -> l2
 rbind(l1, l2) -> l
 l[, offset := cumsum(c(0, length[1:(.N-1)]))]
 l[, plot_offset := cumsum(c(0, length[1:(.N-1)]+1e8))]
 chrNames(agp=T, species)[l, on="agp_chr"] -> l

 copy(assembly_v2$info) -> ai
 ai[!s, on="scaffold"] -> u
 ai[s, on="scaffold"][, scaffold := paste0(scaffold, "_hap1")] -> i1
 ai[s, on="scaffold"][, scaffold := paste0(scaffold, "_hap2")] -> i2
 rbind(u, i1, i2) -> i

 assembly_v2$fpairs[, .(scaffold1, scaffold2, pos1, pos2)] -> f
 f[scaffold1 %in% s, scaffold1 := paste0(scaffold1, ifelse(runif(.N) > 0.5, "_hap1", "_hap2"))]
 f[scaffold2 %in% s, scaffold2 := paste0(scaffold2, ifelse(runif(.N) > 0.5, "_hap1", "_hap2"))]

 list(info=i, fpairs=f) -> assembly_hap
 list(agp=a, chrlen=l) -> hic_map
 list(assembly_hap=assembly_hap, hic_map=hic_map)
}

##
# Assign contigs that are not place in the guide map to haplotype using Hi-C
##

#read haplotype separation and guide map table

readRDS('FB20_005_1_assembly_v2.Rds') -> assembly_v2
readRDS('FB20_005_1_assembly_v2_cov.Rds') -> cov

##
# Construct first Hi-C map
##

#read Hi-C fragment data
f <- '../FB20_005_1_hifiasm.p_utg_MboI_fragments_30bp.bed'
read_fragdata(info=assembly_v2$info, file=f)->frag_data

#read haplotype information
readRDS(file="FB20_005_1_assembly_v2_Hv_guide+HiClift+hap.Rds") -> hh

#Hi-C map for haplotype 1
hh[hap %in% c(1,3), .(scaffold, length, chr=as.integer(substr(chr, 4, 4)), cM = pos/1e6)] -> hic_info
frag_data$info[, .(scaffold, nfrag)][hic_info, on="scaffold"] -> hic_info

hic_map(info=hic_info, assembly=assembly_v2, frags=frag_data$bed, species="barley", ncores=21,
  min_nfrag_scaffold=30, max_cM_dist = 1000,
  binsize=2e5, min_nfrag_bin=10, gap_size=100)->hic_map_v1_hap1

saveRDS(hic_map_v1_hap1, file="FB20_005_1_hic_map_v1_hap1.Rds")

#Hi-C map for haplotype 2
hh[hap %in% c(2,3), .(scaffold, length, chr=as.integer(substr(chr, 4, 4)), cM = pos/1e6)] -> hic_info
frag_data$info[, .(scaffold, nfrag)][hic_info, on="scaffold"] -> hic_info

hic_map(info=hic_info, assembly=assembly_v2, frags=frag_data$bed, species="barley", ncores=21,
  min_nfrag_scaffold=30, max_cM_dist = 1000,
  binsize=2e5, min_nfrag_bin=10, gap_size=100)->hic_map_v1_hap2

saveRDS(hic_map_v1_hap2, file="FB20_005_1_hic_map_v1_hap2.Rds")

#create Hi-C plots for haplotype 1
snuc <- '/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/FB20_005_1/FB20_005_1_hifiasm.p_utg_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="FB20_005_1_hic_map_v1_hap1.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v1_hap1_l

#create Hi-C plots for haplotype 2
snuc <- '/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/FB20_005_1/FB20_005_1_hifiasm.p_utg_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="FB20_005_1_hic_map_v1_hap2.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v1_hap2_l

#combine both Hi-C maps
combine_hic(hic_map_v1_hap1, hic_map_v1_hap2, assembly_v2) -> ch
ch$hic_map -> hic_map_v1
ch$assembly -> assembly_v2_hap

#create inter-chromosomal plot for combined hap1+hap2 Hi-C haplotypes
nuc <- '/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/FB20_005_1/FB20_005_1_hifiasm.p_utg_MboI_fragments_30bp_split.nuc.txt'
add_psmol_fpairs(assembly=assembly_v2_hap, hic_map=hic_map_v1, nucfile=nuc)->hic_map_v1_l

bin_hic_step(hic=hic_map_v1_l$links, frags=hic_map_v1_l$frags, binsize=1e6, chrlen=hic_map_v1_l$chrlen, cores=14)->hic_map_v1_l$hic_1Mb

f <- "FB20_005_1_hic_map_v1_interchromosomal.png"
interchromosomal_matrix_plot(hic_map=hic_map_v1_l, file=f, species="hordeum_bulbosum")

#write Hi-C map for editing
write_hic_map(rds="FB20_005_1_hic_map_v1_hap1.Rds", file="FB20_005_1_hic_map_v1_hap1.xlsx", species="barley")
write_hic_map(rds="FB20_005_1_hic_map_v1_hap2.Rds", file="FB20_005_1_hic_map_v1_hap2.xlsx", species="barley")

###### write fasta

readRDS('FB20_005_1_assembly_v2.Rds') -> assembly_v2
readRDS('FB20_005_1_hic_map_v1_hap1.Rds') -> hic_map_v1_hap1

readRDS('FB20_005_1_hic_map_v1_hap2.Rds') -> hic_map_v1_hap2


fasta <- '/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/FB20_005_1/FB20_005_1_hifiasm.p_utg.fasta'
sink("FB20_005_1_pseudomolecules_v1_hap1.log")
compile_psmol(fasta=fasta, output="FB20_005_1_pseudomolecules_v1_hap1", hic_map=hic_map_v1_hap1, assembly=assembly_v2, cores=30)
sink()

sink("FB20_005_1_pseudomolecules_v1_hap2.log")
compile_psmol(fasta=fasta, output="FB20_005_1_pseudomolecules_v1_hap2", hic_map=hic_map_v1_hap2, assembly=assembly_v2, cores=30)
sink()


###Do mannel check for 2 haplotypes

#import Hi-C maps after editing, v1 -> v2
read_hic_map(rds="FB20_005_1_hic_map_v1_hap1.Rds", file="FB20_005_1_hic_map_v1_hap1_edit.xlsx") -> nmap
hic_map(species="barley", agp_only=T, map=nmap)->hic_map_v2_hap1
saveRDS(hic_map_v2_hap1, file="FB20_005_1_hic_map_v2_hap1.Rds")

read_hic_map(rds="FB20_005_1_hic_map_v1_hap2.Rds", file="FB20_005_1_hic_map_v1_hap2_edit.xlsx") -> nmap
hic_map(species="barley", agp_only=T, map=nmap)->hic_map_v2_hap2
saveRDS(hic_map_v2_hap2, file="FB20_005_1_hic_map_v2_hap2.Rds")

#proceed with Hi-C plots and repeat edit cycle if need be


#create Hi-C plots for haplotype 1
snuc <- '/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/FB20_005_1/FB20_005_1_hifiasm.p_utg_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="FB20_005_1_hic_map_v2_hap1.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v2_hap1_l

#create Hi-C plots for haplotype 2
snuc <- '/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/FB20_005_1/FB20_005_1_hifiasm.p_utg_MboI_fragments_30bp_split.nuc.txt'
hic_plots(rds="FB20_005_1_hic_map_v2_hap2.Rds", cov=F, assembly=assembly_v2, cores=30, species="barley", nuc=snuc) -> hic_map_v1_hap2_l



#readRDS('FB20_005_1_assembly_v2.Rds') -> assembly_v2
readRDS('FB20_005_1_hic_map_v2_hap1.Rds') -> hic_map_v2_hap1
readRDS('FB20_005_1_hic_map_v2_hap2.Rds') -> hic_map_v2_hap2


fasta <- '/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/FB20_005_1/FB20_005_1_hifiasm.p_utg.fasta'
sink("FB20_005_1_pseudomolecules_v2_hap1.log")
compile_psmol(fasta=fasta, output="FB20_005_1_pseudomolecules_v2_hap1", hic_map=hic_map_v2_hap1, assembly=assembly_v2, cores=30)
sink()

sink("FB20_005_1_pseudomolecules_v2_hap2.log")
compile_psmol(fasta=fasta, output="FB20_005_1_pseudomolecules_v2_hap2", hic_map=hic_map_v2_hap2, assembly=assembly_v2, cores=30)
sink()


combine_hic(hic_map_v2_hap1, hic_map_v2_hap2, assembly_v2) -> ch2
ch2$hic_map -> hic_map_v2
ch2$assembly -> assembly_v2_hap

nuc <- '/filer-dg/agruppen/dg7/fengj/hbulbosum_genome/FB20_005_1/FB20_005_1_hifiasm.p_utg_MboI_fragments_30bp_split.nuc.txt'
add_psmol_fpairs(assembly=assembly_v2_hap, hic_map=hic_map_v2, nucfile=nuc)->hic_map_v2_l

bin_hic_step(hic=hic_map_v2_l$links, frags=hic_map_v2_l$frags, binsize=1e6, chrlen=hic_map_v2_l$chrlen, cores=14)->hic_map_v2_l$hic_1Mb

f <- "FB20_005_1_hic_map_v2_interchromosomal.png"
interchromosomal_matrix_plot(hic_map=hic_map_v2_l, file=f, species="hordeum_bulbosum")

