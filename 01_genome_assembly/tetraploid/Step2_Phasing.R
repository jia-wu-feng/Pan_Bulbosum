#https://bitbucket.org/tritexassembly/tritexassembly.bitbucket.io/src/master/R/
source('/filer-dg/agruppen/seq_shared/mascher/code_repositories/triticeae.bitbucket.io/R/pseudomolecule_construction.R')

#read chromosome lengths of Morex V3
fread('/filer-dg/agruppen/seq_shared/mascher/morexV3_pseudomolecules_200421/200416_MorexV3_pseudomolecules.fasta.fai', sel=1:2, col.names=c("chr", "len"))->fai

#read centromere positions
fread('/filer-dg/agruppen/seq_shared/mascher/morexV3_pseudomolecules_200421/MorexV3_centromere_positions.tsv')->cen
setnames(cen, c("chr", "cen_pos"))

#read Hi-C mapping output (list of read pairs connecting RE fragments)
fread('/filer-dg/agruppen/dg2/fengj/hbulbosum_genome/genome_assembly/A40/hic/hicreads/hicreads_fragment_pairs.tsv.gz')->f
setnames(f, c("ctg1", "pos1", "ctg2", "pos2"))

#read contig lengths of assebmbly
faif <- '/filer-dg/agruppen/dg2/fengj/hbulbosum_genome/genome_assembly/A40/A40_hifiasm.fasta.fai'
fread(faif, sel=1:2, col.names=c("contig", "contig_length"))->fb

#read assembly object. This is the hifiasm unitig assembly after breaking chimeras based on coverage drops
readRDS('A40_assembly_v2.Rds') -> assembly_v2

##
# Positions contigs based on guide map
##

#read alignment of MorexV3 HC genes 
fread('./A40_hifiasm_MorexV3_HC_table.txt')->x
setnames(x, c("contig", "start", "end", "transcript", "alnlen", "id"))

#read positions of MorexV3 genes
fread('/filer-dg/agruppen/seq_shared/mascher/morex_v3_annotation_200818/Hv_Morex.pgsb.Jul2020.HC_mrna_pos.txt')->p
setnames(p, c("chr", "chr_start", "chr_end", "transcript"))

#keep only genes aligned for >= 90 % of their length with 90 % identity, allow up to two alignment (one per haplotype), modify for tetraploids
p[x, on="transcript"] -> px
px[alnlen >= 90 & id >= 90]->a
px[a[, .N, key=transcript][N <= 4]$transcript, on="transcript"]->aa

#convert original contig positions to corrected (assembly_v2) positions 
assembly_v2$info[, .(contig=orig_scaffold, start=orig_start, orig_start, scaffold)][aa, on=c("contig", "start"), roll=T] -> aa

#chromosome assignment
aa[, .N, key=.(chr, scaffold)][, p := N/sum(N), by=scaffold][order(-p)][!duplicated(scaffold)] -> cc

#keep contigs with at least 3 aligned genes, 65 % of aligned are from the major chromosome
cc[N >= 3 & p >= 0.65] -> cc

#check how much of the assembly can be assigned to chromosomes
assembly_v2$info[, .(scaffold, scaffold_length=length)][cc, on="scaffold"] -> cc
cc[scaffold_length>=1e6]->cc 
sum(cc$scaffold_length)/1e6

####11865.57 


#get the approximate chromosome positions (median of alignment coordinates)
aa[cc[, .(scaffold, chr)], on=c("chr", "scaffold")][, .(pos=median(as.numeric(chr_start)), pos_mad=mad(as.numeric(chr_start))), key=scaffold][cc, on="scaffold"] -> cc
cc[, mr := pos_mad / scaffold_length]  
setorder(cc, chr, pos) 
cc[, chr_idx := 1:.N, by=chr]
cc[, agp_pos := c(0, cumsum(scaffold_length[-.N])), by=chr]

#save results
saveRDS(cc, file="A40_assembly_v2_Hv_guide.Rds")

##
# Positioning additional contigs by Hi-C
## 

#load guide map positions
readRDS('A40_assembly_v2_Hv_guide.Rds') -> cc
#get Hi-C links in the terminal 2 Mb of each scaffold
assembly_v2$fpairs[, .(scaffold=scaffold1, pos=pos1, link=scaffold2)] -> ff
assembly_v2$info[, .(scaffold, length)][ff, on="scaffold"] -> ff
ff[pos <= 2e6 | length - pos <= 2e6] -> ff
ff[scaffold != link] -> ff
ff[, .N, key=.(scaffold, link)] -> fa
#exclude scaffold-link paris with only a single Hi-C pair
fa[N > 1] -> fa
#add guide map positions
cc[, .(scaffold, scaffold_chr=chr, scaffold_pos=pos)][fa, on="scaffold"] -> fa
cc[, .(link=scaffold, link_chr=chr, link_pos=pos)][fa, on="link"] -> fa
#get approximate position based on Hi-C links
fa[!is.na(link_chr), .(n=sum(N), pos=weighted.mean(link_pos, N)), key=.(scaffold, scaffold_chr, scaffold_pos, link_chr)] -> fv
fv[, p := n/sum(n), by=scaffold]
assembly_v2$info[, .(scaffold, length)][fv, on="scaffold"] -> fv
#get chromosome assignment for scaffolds without position in guide map
fv[is.na(scaffold_chr)][order(-p)][!duplicated(scaffold)][order(-length)] -> cc_lift0

#exclude contigs shorter than 300 kb
cc_lift0[length >= 3e5] -> cc_lift
#merge guide map and Hi-C lift tables; write output
rbind(cc[, .(scaffold, chr, pos)], cc_lift[, .(scaffold, chr=link_chr, pos)]) -> a
a[assembly_v2$info[, .(scaffold, length)], on="scaffold"] -> a
saveRDS(a, file="A40_assembly_v2_Hv_guide+HiClift.Rds")

#check how much of the assembly is placed
 sum(a[!is.na(chr)]$length)
#[1] 1.348e+10
 sum(a$length)
#[1] 14236991145
 1.348e+10/14236991145
#[1] 0.9468293

##
# Haplotype phasing
##

#read lengths of Morex chromosomes
fread('/filer-dg/agruppen/seq_shared/mascher/morexV3_pseudomolecules_200421/200416_MorexV3_pseudomolecules.fasta.fai', sel=1:2, col.names=c("chr", "len"))->morexfai

#read centromere position in Morex 
fread('/filer-dg/agruppen/seq_shared/mascher/morexV3_pseudomolecules_200421/MorexV3_centromere_positions.tsv')->cen
setnames(cen, c("chr", "cen_pos"))

#read scaffold posiitions and add them to the Hi-C link table
readRDS(file="A40_assembly_v2_Hv_guide.Rds") -> cc
assembly_v2$fpairs[, .N, key=.(ctg1=scaffold1, ctg2=scaffold2)] -> v
cc[, .(ctg1=scaffold, chr1=chr, pos1=pos)][v, on="ctg1"] -> v
cc[, .(ctg2=scaffold, chr2=chr, pos2=pos)][v, on="ctg2"] -> v

#run a PCA on the intra-chromosomal matrices
rbindlist(mclapply(mc.cores=7, morexfai[1:7, chr], function(j){
 #fill in empty scaffold pairs with 0
 setnames(v[chr1 == chr2 & chr1 == j][, chr2 := NULL], "chr1", "chr") -> x
 dcast(x, ctg1 ~ ctg2, value.var="N", fill = 0) -> x
 melt(x, id="ctg1", measure=setdiff(names(x), "ctg1"), variable.factor=F, value.name="N", variable.name="ctg2") -> x
 x[, l := log10(0.1 + N)]
 cc[, .(ctg1=scaffold, chr, pos1=pos)][x, on="ctg1"] -> x
 cc[, .(ctg2=scaffold, pos2=pos)][x, on="ctg2"] -> x
 morexfai[x, on="chr"] -> x
 cen[x, on="chr"] -> x
 x[, end_dist1 := pmin(pos1, len - pos1)]
 x[, end_dist2 := pmin(pos2, len - pos2)]
 x[, cen_dist1 := abs(cen_pos - pos1)]
 x[, cen_dist2 := abs(cen_pos - pos2)]
 x[, ldist := abs(pos1 - pos2)]
 #get "normalized" Hi-C ounts by removing the factor linear distance between loci, distance from centromere and distance from chromosome end (all log scaled)
 x[, res := lm(l ~ log(1+ldist) * log(cen_dist1) * log(cen_dist2) * log(end_dist1) * log(end_dist2))$res]

 #convert to matrix
 dcast(x, ctg1 ~ ctg2, value.var="res", fill=0) -> y
 y[, ctg1 := NULL]
 #run PCA on correlation matrix
 prcomp(cor(y), scale=T, center=T)->pca
 #get first four eigenvector
 data.table(contig=rownames(pca$rotation), pca$rotation[, 1:4]) -> p
 setorder(cc[, .(contig=scaffold, chr, pos)][p, on="contig"], chr, pos) -> pp
 pp[, chr := j]
})) -> pp

assembly_v2$info[, .(contig=scaffold, contig_length=length)][pp, on="contig"] -> pp
pp[, idx := 1:.N]
#save results
saveRDS(pp, file="A40_assembly_v2_haplotype_separation_HiC.Rds")

#read coverage information and convert to correct assembly coordinates
fread(cmd="grep '^S' ../A40_hifiasm.p_utg.noseq.gfa | cut -f 2,5 | tr ':' '\\t' | cut -f 1,4", head=F) -> cov
setnames(cov, c("contig", "cc"))
cov[, contig := sub('l$', '', sub("utg0*", "contig_", contig))]
assembly_v2$info[, .(contig=orig_scaffold, scaffold)][cov, on=c("contig")] -> cov
setorder(cc[, .(scaffold, scaffold_length, chr, pos)][cov, on="scaffold"], chr, pos) -> cov
cov[, idx := 1:.N]
saveRDS(file="A40_assembly_v2_cov.Rds", cov)


cov[, .(contig=scaffold, cc)][pp, on="contig"] -> pp

pp[cc<45, col := "red"]
pp[cc>=45&cc<80,col:="blue"]
pp[cc>=80&cc<120,col:="black"]
pp[cc>=120,col:="orange"]


pdf("A40_haplotype_separation_HiC_all_pc1_2.pdf", height=10, width=10)

lapply(morexfai[1:7, chr], function(i){
  pp[chr == i, plot(PC2, PC1, col=col, pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC1",
            xlab="PC2")]
  legend("topleft", pch=19, bty='n', legend=c("1X", "2X","3X","4X"), col=c("red", "blue","black","orange"))
 })

dev.off()

pdf("A40_haplotype_separation_HiC_all_pc2_3.pdf", height=10, width=10)

lapply(morexfai[1:7, chr], function(i){
  pp[chr == i, plot(PC3, PC2, col=col, pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC2",
            xlab="PC3")]
  legend("topleft", pch=19, bty='n', legend=c("1X", "2X","3X","4X"), col=c("red", "blue","black","orange"))
 })

dev.off()

pdf("A40_haplotype_separation_HiC_all_pc3_4.pdf", height=10, width=10)

lapply(morexfai[1:7, chr], function(i){
  pp[chr == i, plot(PC4, PC3, col=col, pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC3",
            xlab="PC4")]
  legend("topleft", pch=19, bty='n', legend=c("1X", "2X","3X","4X"), col=c("red", "blue","black","orange"))
 })

dev.off()


pdf("A40_haplotype_separation_HiC_all.pdf", height=4, width=10)
par(mfrow=c(2,4))
lapply(c("PC1", "PC2", "PC3", "PC4","cc"), function(j){
 lapply(morexfai[1:7, chr], function(i){
  pp[chr == i, plot(pos/1e6, las=1, bty='l', get(j), type='n', main=sub("chr", "", chr[1]), ylab=j,
            xlab="Hv syntenic position [Mb]")]
  pp[chr == i, lines(lwd=3, c(pos/1e6, (pos + contig_length)/1e6), c(get(j), get(j))), by=idx]
 })
 plot(axes=F, xlab="", ylab="", 0, type='n')
})
dev.off()




### filter coverage 45 assembly with 2X coverage

rbindlist(mclapply(mc.cores=7, morexfai[1:7, chr], function(j){
 #fill empty values in link list (with 0)
 setnames(v[chr1 == chr2 & chr1 == j][, chr2 := NULL], "chr1", "chr") -> x
 x[ctg1 %in%  cov[cc <= 45]$scaffold & ctg2 %in% cov[cc <= 45]$scaffold] -> x
 dcast(x, ctg1 ~ ctg2, value.var="N", fill = 0) -> x
 melt(x, id="ctg1", measure=setdiff(names(x), "ctg1"), variable.factor=F, value.name="N", variable.name="ctg2") -> x
 #convert link counts to log
 x[, l := log10(0.1 + N)]
 #merge with positional information 
 cc[, .(ctg1=scaffold, chr, pos1=pos)][x, on="ctg1"] -> x
 cc[, .(ctg2=scaffold, pos2=pos)][x, on="ctg2"] -> x
 morexfai[x, on="chr"] -> x
 cen[x, on="chr"] -> x
 #specify linear model to remove DDD and Rabl effects 
 x[, end_dist1 := pmin(pos1, len - pos1)]
 x[, end_dist2 := pmin(pos2, len - pos2)]
 x[, cen_dist1 := abs(cen_pos - pos1)]
 x[, cen_dist2 := abs(cen_pos - pos2)]
 x[, ldist := abs(pos1 - pos2)]
 x[, res := lm(l ~ log(1+ldist) * log(cen_dist1) * log(cen_dist2) * log(end_dist1) * log(end_dist2))$res]

 #run PCA on residuals
 dcast(x, ctg1 ~ ctg2, value.var="res", fill=0) -> y
 y[, ctg1 := NULL]
 #PCA of correlation matrix, e.g. SVD of co-variance of correlation of Hi-C distance matrix (seems complicated, but gives good results)
 prcomp(cor(y), scale=T, center=T)->pca
 #extract first four PCs
 data.table(contig=rownames(pca$rotation), pca$rotation[, 1:4]) -> p
 #merge PC table with positional information
 setorder(cc[, .(contig=scaffold, chr, pos)][p, on="contig"], chr, pos) -> pp4
 pp4[, chr := j]
})) -> pp4

assembly_v2$info[, .(contig=scaffold, contig_length=length)][pp4, on="contig"] -> pp4

pp4[, idx := 1:.N]


cov[, .(contig=scaffold, cc)][pp4, on="contig"] -> pp4

pp4[cc<45, col := "red"]
pp4[cc>=45&cc<80,col:="blue"]
pp4[cc>=80&cc<120,col:="black"]
pp4[cc>=120,col:="orange"]


pdf("A40_haplotype_separation_HiC_filter_coverage_pc1_2.pdf", height=10, width=10)

lapply(morexfai[1:7, chr], function(i){
  pp4[chr == i, plot(PC2, PC1, col=col, pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC1",
            xlab="PC2")]
   legend("topleft", pch=19, bty='n', legend=c("1X", "2X","3X","4X"), col=c("red", "blue","black","orange"))
 })

dev.off()

pdf("A40_haplotype_separation_HiC_filter_coverage_pc2_3.pdf", height=10, width=10)

lapply(morexfai[1:7, chr], function(i){
  pp4[chr == i, plot(PC3, PC2, col=col, pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC2",
            xlab="PC3")]
  legend("topleft", pch=19, bty='n', legend=c("1X", "2X","3X","4X"), col=c("red", "blue","black","orange"))
 })

dev.off()

pdf("A40_haplotype_separation_HiC_filter_coverage_pc3_4.pdf", height=10, width=10)

lapply(morexfai[1:7, chr], function(i){
  pp4[chr == i, plot(PC4, PC3, col=col, pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC3",
            xlab="PC4")]
  legend("topleft", pch=19, bty='n', legend=c("1X", "2X","3X","4X"), col=c("red", "blue","black","orange"))
 })

dev.off()


saveRDS(pp4, file="A40_assembly_v2_haplotype_separation_HiC_all.Rds")






#Pre-cluster PCA results using Kmeans

readRDS(file="A40_assembly_v2_haplotype_separation_HiC_all.Rds") -> pp4

pp4[,cluster:=0]
primer_cluster <- function(pp4, chrj, center){
data.frame(pp4[chr==chrj,.(PC1,PC2,PC3,PC4)])-> mm
kmeans(mm,centers=center, iter.max = 30, nstart = 30)$cluster->q1
data.table(pp4[chr==chrj])-> mm1
mm1[,cluster:=q1]
rbind(pp4[chr!=chrj], mm1)->pp4
return(pp4) 
}

primer_cluster(pp4,"chr1H",4)-> pp4
primer_cluster(pp4,"chr2H",4)-> pp4
primer_cluster(pp4,"chr3H",4)-> pp4
primer_cluster(pp4,"chr4H",4)-> pp4
primer_cluster(pp4,"chr5H",4)-> pp4
primer_cluster(pp4,"chr6H",4)-> pp4
primer_cluster(pp4,"chr7H",4)-> pp4

pp[,cluster:=0]->pp1
rbind( pp1[!pp4$contig, on="contig"][, PC1 := 0 ][, PC2 := 0 ][, PC3 := 0 ][, PC4 := 0 ], pp4)->pq

### Distinguish mixed haplotypes using coverage

pq[contig=="contig_corrected_v1_5136",cluster := 22]
pq[contig=="contig_corrected_v1_95",cluster := 22]
pq[cc>=45 & cc<80, cluster := 22]
pq[cc>=80 & cc<120, cluster := 33]
pq[cc>=120, cluster :=44]

pq[, col := "black"]
pq[cluster == 1, col := "red"]
pq[cluster == 2, col := "blue"]
pq[cluster == 3, col := "purple"]
pq[cluster == 4, col := "orange"]

pq[cluster == 22 |cluster == 33 | cluster == 44, col :="black"]


setorder(pq, chr, pos)->pq
pq[, idx := 1:.N]

saveRDS(pq, file="A40_assembly_v2_manual.Rds")



pdf("A40_haplotype_separation_HiC_filter_coverage_45_pc1_2_cluster.pdf", height=10, width=10)

lapply(morexfai[1:7, chr], function(i){
  pq[chr == i, plot(PC2, PC1, col=cluster, pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC1",
            xlab="PC2")]
  legend("topleft", pch=19, bty='n', legend=c("1", "2","3","4"), col=c(1,2,3,4))
 })

dev.off()

pdf("A40_haplotype_separation_HiC_filter_coverage_45_pc2_3_cluster.pdf", height=10, width=10)

lapply(morexfai[1:7, chr], function(i){
  pq[chr == i, plot(PC3, PC2, col=cluster, pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC2",
            xlab="PC3")]
  legend("topleft", pch=19, bty='n', legend=c("1", "2","3","4"), col=c(1,2,3,4))
 })

dev.off()

pdf("A40_haplotype_separation_HiC_filter_coverage_45_pc3_4_cluster.pdf", height=10, width=10)

lapply(morexfai[1:7, chr], function(i){
  pq[chr == i, plot(PC4, PC3, col=cluster, pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC3",
            xlab="PC4")]
  legend("topleft", pch=19, bty='n', legend=c("1", "2","3","4"), col=c(1,2,3,4))
 })

dev.off()

### Calculate the length of the cluster based on coverage

cal_chrlen <- function(pq, chrj,cov1){
  hap1=sum(pq[chr == chrj & cluster==1 & cc<=cov1]$contig_length)+sum(pq[chr == chrj & cluster==1 & cc>=cov1]$contig_length)*2
  hap2=sum(pq[chr == chrj & cluster==2 & cc<=cov1]$contig_length)+sum(pq[chr == chrj & cluster==2 & cc>=cov1]$contig_length)*2
  hap3=sum(pq[chr == chrj & cluster==3 & cc<=cov1]$contig_length)+sum(pq[chr == chrj & cluster==3 & cc>=cov1]$contig_length)*2
  hap4=sum(pq[chr == chrj & cluster==4 & cc<=cov1]$contig_length)+sum(pq[chr == chrj & cluster==4 & cc>=cov1]$contig_length)*2
  print(hap1)
  print(hap2)
  print(hap3)
  print(hap4)
}

cal_chrlen(pq,"chr1H",45)
cal_chrlen(pq,"chr2H",45)
cal_chrlen(pq,"chr3H",45)
cal_chrlen(pq,"chr4H",45)
cal_chrlen(pq,"chr5H",45)
cal_chrlen(pq,"chr6H",45)
cal_chrlen(pq,"chr7H",45)

> cal_chrlen(pq,"chr1H",45)
[1] 394050661
[1] 354501680
[1] 379430333
[1] 347908367
> cal_chrlen(pq,"chr2H",45)
[1] 449199099
[1] 451085238
[1] 504001382
[1] 461401651
> cal_chrlen(pq,"chr3H",45)
[1] 456338782
[1] 474787619
[1] 478076326
[1] 483633164
> cal_chrlen(pq,"chr4H",45)
[1] 381461783
[1] 417103232
[1] 397240082
[1] 388591793
> cal_chrlen(pq,"chr5H",45)
[1] 447148528
[1] 412256251
[1] 414050167
[1] 418149355
> cal_chrlen(pq,"chr6H",45)
[1] 374488285
[1] 345488812
[1] 336125826
[1] 370503325
> cal_chrlen(pq,"chr7H",45)
[1] 459964897
[1] 480384645
[1] 440080902
[1] 423418976


###cluster Re- PCA function 

###Repca is used to distinguish clusters more than 1

Repca <- function(pq, chrj, cov,v, cl, plevel=0)
{
pq=pq
j=chrj
setnames(v[chr1 == chr2 & chr1 == j][, chr2 := NULL], "chr1", "chr") -> x
x[ctg1 %in%  pq[cluster==cl ]$contig & ctg2 %in% pq[cluster==cl ]$contig] -> x
dcast(x, ctg1 ~ ctg2, value.var="N", fill = 0) -> x
melt(x, id="ctg1", measure=setdiff(names(x), "ctg1"), variable.factor=F, value.name="N", variable.name="ctg2") -> x
x[, l := log10(0.1 + N)]
cc[, .(ctg1=scaffold, chr, pos1=pos)][x, on="ctg1"] -> x
cc[, .(ctg2=scaffold, pos2=pos)][x, on="ctg2"] -> x
morexfai[x, on="chr"] -> x
cen[x, on="chr"] -> x
x[, end_dist1 := pmin(pos1, len - pos1)]
x[, end_dist2 := pmin(pos2, len - pos2)]
x[, cen_dist1 := abs(cen_pos - pos1)]
x[, cen_dist2 := abs(cen_pos - pos2)]
x[, ldist := abs(pos1 - pos2)]
x[, res := lm(l ~ log(1+ldist) * log(cen_dist1) * log(cen_dist2) * log(end_dist1) * log(end_dist2))$res]
dcast(x, ctg1 ~ ctg2, value.var="res", fill=0) -> y
y[, ctg1 := NULL]
prcomp(cor(y), scale=T, center=T)->pca
data.table(contig=rownames(pca$rotation), pca$rotation[, 1:4]) -> p
setorder(cc[, .(contig=scaffold, chr, pos)][p, on="contig"], chr, pos) -> pp1c1
pp1c1[, chr := j]
assembly_v2$info[, .(contig=scaffold, contig_length=length)][pp1c1, on="contig"] -> pp1c1
pp1c1[, idx := 1:.N]
cov[, .(contig=scaffold, cc)][pp1c1, on="contig"] -> pp1c1
pp1c1[cc<32, col := "red"]
pp1c1[cc>=32&cc<64,col:="blue"]
pp1c1[cc>=64&cc<96,col:="black"]
pp1c1[cc>=96,col:="orange"]
cov[, .(contig=scaffold, cc)][pp1c1, on="contig"] -> pp1c1
pp1c1[PC1>plevel, cluster := cl]
pp1c1[PC1<plevel, cluster := 4]
pdf(paste("haplotype_separation_HiC_pp",chrj,"c",cl,".pdf",sep=""), height=8, width=8)
plot(pp1c1$PC2, pp1c1$PC1, col=pp1c1$col, pch=20, las=1,bty='l', type="p", main=paste(chrj," cluster",cl,sep=""), ylab="PC1",xlab="PC2")
dev.off()
rbind(pq[!pp1c1$contig, on="contig"], pp1c1)-> pq1
print(sum(pq1[chr == chrj & cluster==1]$contig_length))
print(sum(pq1[chr == chrj & cluster==2]$contig_length))
print(sum(pq1[chr == chrj & cluster==3]$contig_length))
print(sum(pq1[chr == chrj & cluster==4]$contig_length))
pp1c1[, idx := 1:.N]
pdf(paste("haplotype_separation_HiC_pp",chrj,"c",cl,"_repca.pdf",sep=""), height=4, width=10)
lapply(c("PC1", "PC2", "PC3", "PC4"), function(j){
 pp1c1[chr==chrj, plot(pos/1e6, las=1, bty='l', get(j), type='n', main=sub("chr", "", chrj), ylab=j,
            xlab="Hv syntenic position [Mb]")]
 pp1c1[chr==chrj, lines(lwd=3, c(pos/1e6, (pos + contig_length)/1e6), c(get(j), get(j))), by=idx]
})
dev.off()
return(pq1) 
}

### Repcafor2 is designed to redistinguish and demarcate the boundaries between two haplotypes.

Repcafor2 <- function(pq, chrj, cov,v, cl,cl2, plevel=0)
{
pq=pq
j=chrj
setnames(v[chr1 == chr2 & chr1 == j][, chr2 := NULL], "chr1", "chr") -> x
x[ctg1 %in%  pq[cluster==cl | cluster==cl2 ]$contig & ctg2 %in% pq[cluster==cl | cluster==cl2 ]$contig] -> x
dcast(x, ctg1 ~ ctg2, value.var="N", fill = 0) -> x
melt(x, id="ctg1", measure=setdiff(names(x), "ctg1"), variable.factor=F, value.name="N", variable.name="ctg2") -> x
x[, l := log10(0.1 + N)]
cc[, .(ctg1=scaffold, chr, pos1=pos)][x, on="ctg1"] -> x
cc[, .(ctg2=scaffold, pos2=pos)][x, on="ctg2"] -> x
morexfai[x, on="chr"] -> x
cen[x, on="chr"] -> x
x[, end_dist1 := pmin(pos1, len - pos1)]
x[, end_dist2 := pmin(pos2, len - pos2)]
x[, cen_dist1 := abs(cen_pos - pos1)]
x[, cen_dist2 := abs(cen_pos - pos2)]
x[, ldist := abs(pos1 - pos2)]
x[, res := lm(l ~ log(1+ldist) * log(cen_dist1) * log(cen_dist2) * log(end_dist1) * log(end_dist2))$res]
dcast(x, ctg1 ~ ctg2, value.var="res", fill=0) -> y
y[, ctg1 := NULL]
prcomp(cor(y), scale=T, center=T)->pca
data.table(contig=rownames(pca$rotation), pca$rotation[, 1:4]) -> p
setorder(cc[, .(contig=scaffold, chr, pos)][p, on="contig"], chr, pos) -> pp1c1
pp1c1[, chr := j]
assembly_v2$info[, .(contig=scaffold, contig_length=length)][pp1c1, on="contig"] -> pp1c1
pp1c1[, idx := 1:.N]
cov[, .(contig=scaffold, cc)][pp1c1, on="contig"] -> pp1c1
pp1c1[cc<32, col := "red"]
pp1c1[cc>=32&cc<64,col:="blue"]
pp1c1[cc>=64&cc<96,col:="black"]
pp1c1[cc>=96,col:="orange"]
cov[, .(contig=scaffold, cc)][pp1c1, on="contig"] -> pp1c1
pp1c1[, cluster := 5]
pdf(paste("haplotype_separation_HiC_pp",chrj,"c",cl,cl2,".pdf",sep=""), height=8, width=8)
plot(pp1c1$PC2, pp1c1$PC1, col=pp1c1$col, pch=20, las=1,bty='l', type="p", main=paste(chrj," cluster",cl,sep=""), ylab="PC1",xlab="PC2")
dev.off()
rbind(pq[!pp1c1$contig, on="contig"], pp1c1)-> pq1
print(sum(pq1[chr == chrj & cluster==1]$contig_length))
print(sum(pq1[chr == chrj & cluster==2]$contig_length))
print(sum(pq1[chr == chrj & cluster==3]$contig_length))
print(sum(pq1[chr == chrj & cluster==4]$contig_length))
print(sum(pq1[chr == chrj & cluster==5]$contig_length))
pp1c1[, idx := 1:.N]
pdf(paste("haplotype_separation_HiC_pp",chrj,"c",cl,cl2,"_repca.pdf",sep=""), height=4, width=10)
lapply(c("PC1", "PC2", "PC3", "PC4"), function(j){
 pp1c1[chr==chrj, plot(pos/1e6, las=1, bty='l', get(j), type='n', main=sub("chr", "", chrj), ylab=j,
            xlab="Hv syntenic position [Mb]")]
 pp1c1[chr==chrj, lines(lwd=3, c(pos/1e6, (pos + contig_length)/1e6), c(get(j), get(j))), by=idx]
})
dev.off()
return(pq1) 
}

### Repcafor2check is used to check whether there is a clear separation between two haplotypes. 

Repcafor2check <- function(pq, chrj, cov,v, cl,cl2, prefix, plevel=0,pc=1,error_length=0)
{
pq=pq
j=chrj
dir.create(prefix)
lencl=sum(pq[chr == chrj & cluster==cl]$contig_length)
lencl2=sum(pq[chr == chrj & cluster==cl2]$contig_length)
setnames(v[chr1 == chr2 & chr1 == j][, chr2 := NULL], "chr1", "chr") -> x
x[ctg1 %in%  pq[cluster==cl | cluster==cl2 ]$contig & ctg2 %in% pq[cluster==cl | cluster==cl2 ]$contig] -> x
dcast(x, ctg1 ~ ctg2, value.var="N", fill = 0) -> x
melt(x, id="ctg1", measure=setdiff(names(x), "ctg1"), variable.factor=F, value.name="N", variable.name="ctg2") -> x
x[, l := log10(0.1 + N)]
cc[, .(ctg1=scaffold, chr, pos1=pos)][x, on="ctg1"] -> x
cc[, .(ctg2=scaffold, pos2=pos)][x, on="ctg2"] -> x
morexfai[x, on="chr"] -> x
cen[x, on="chr"] -> x
x[, end_dist1 := pmin(pos1, len - pos1)]
x[, end_dist2 := pmin(pos2, len - pos2)]
x[, cen_dist1 := abs(cen_pos - pos1)]
x[, cen_dist2 := abs(cen_pos - pos2)]
x[, ldist := abs(pos1 - pos2)]
x[, res := lm(l ~ log(1+ldist) * log(cen_dist1) * log(cen_dist2) * log(end_dist1) * log(end_dist2))$res]
dcast(x, ctg1 ~ ctg2, value.var="res", fill=0) -> y
y[, ctg1 := NULL]
prcomp(cor(y), scale=T, center=T)->pca
data.table(contig=rownames(pca$rotation), pca$rotation[, 1:4]) -> p
setorder(cc[, .(contig=scaffold, chr, pos)][p, on="contig"], chr, pos) -> pp1c1
pp1c1[, chr := j]
assembly_v2$info[, .(contig=scaffold, contig_length=length)][pp1c1, on="contig"] -> pp1c1
pp1c1[, idx := 1:.N]
cov[, .(contig=scaffold, cc)][pp1c1, on="contig"] -> pp1c1
pp1c1[cc<32, col := "red"]
pp1c1[cc>=32&cc<64,col:="blue"]
pp1c1[cc>=64&cc<96,col:="black"]
pp1c1[cc>=96,col:="orange"]
cov[, .(contig=scaffold, cc)][pp1c1, on="contig"] -> pp1c1
if(pc==1){
  pp1c1[PC1>plevel, cluster := cl]
  pp1c1[PC1<plevel, cluster := cl2]
}else if(pc==2){
  pp1c1[PC2>plevel, cluster := cl]
  pp1c1[PC2<plevel, cluster := cl2]
}else if(pc==3){
  pp1c1[PC3>plevel, cluster := cl]
  pp1c1[PC3<plevel, cluster := cl2]
}else if(pc==4){
  pp1c1[PC4>plevel, cluster := cl]
  pp1c1[PC4<plevel, cluster := cl2]
}
pdf(paste(prefix,"/haplotype_separation_HiC_pp",chrj,"c",cl,cl2,"_check.pdf",sep=""), height=8, width=8)
plot(pp1c1$PC2, pp1c1$PC1, col=pp1c1$col, pch=20, las=1,bty='l', type="p", main=paste(chrj," cluster",cl,sep=""), ylab="PC1",xlab="PC2")
dev.off()
pp1c1[, idx := 1:.N]
pp1c1->ppp
pq[, .(contig, clo=cluster)][ppp, on="contig"] -> ppp
pdf(paste(prefix,"/haplotype_separation_HiC_pp",chrj,"c",cl,cl2,"_repca_ckeck.pdf",sep=""), height=4, width=10)
lapply(c("PC1", "PC2", "PC3", "PC4"), function(j){
 ppp[chr==chrj, plot(pos/1e6, las=1, bty='l', get(j), type='n', main=sub("chr", "", chrj), ylab=j,
            xlab="Hv syntenic position [Mb]")]
 ppp[chr==chrj, lines(lwd=3, c(pos/1e6, (pos + contig_length)/1e6), c(get(j), get(j)),col=clo), by=idx]
 abline(h=plevel, col="blue")
})
dev.off()
rbind(pq[!pp1c1$contig, on="contig"], pp1c1)-> pq1
print(sum(pq1[chr == chrj & cluster==1]$contig_length))
print(sum(pq1[chr == chrj & cluster==2]$contig_length))
print(sum(pq1[chr == chrj & cluster==3]$contig_length))
print(sum(pq1[chr == chrj & cluster==4]$contig_length))
print(sum(pq1[chr == chrj & cluster==5]$contig_length))
checklencl=sum(pq1[chr == chrj & cluster==cl]$contig_length)
checklencl2=sum(pq1[chr == chrj & cluster==cl2]$contig_length)
if((abs(lencl - checklencl) <= error_length) | (abs(lencl - checklencl2) <= error_length)){
    print("Correct")
    print(abs(lencl - checklencl))
    print(abs(lencl - checklencl2))
}else{
    print("Need re-pca or check plevel")
}
}

####RE-PCA and check for 2

readRDS(file="A40_assembly_v2_manual.Rds") -> pq1

cov[, .(contig=scaffold, cc)][pq1, on="contig"] -> pq1

##check chr1H

# re-assignment

Repcafor2(pq1,"chr1H",cov,v,2,4,0) -> pq1

pq1[chr=="chr1H"&cluster==5& PC1 >= 0 , cluster:=2]
pq1[chr=="chr1H"&cluster==5& PC1 <= 0 , cluster:=4]

Repcafor2check(pq1,"chr1H",cov,v,1,2,"chr1Hcheck",0)

Repcafor2check(pq1,"chr1H",cov,v,1,3,"chr1Hcheck",0)

Repcafor2check(pq1,"chr1H",cov,v,1,4,"chr1Hcheck",0)

Repcafor2check(pq1,"chr1H",cov,v,2,3,"chr1Hcheck",0)

Repcafor2check(pq1,"chr1H",cov,v,2,4,"chr1Hcheck",0)

Repcafor2check(pq1,"chr1H",cov,v,3,4,"chr1Hcheck",0)


cal_chrlen(pq1,"chr1H",45)


####check chr2H

Repcafor2(pq1,"chr2H",cov,v,1,2,0) -> pq1

pq1[chr=="chr2H"&cluster==5& PC1 >= -0.005 , cluster:=1]
pq1[chr=="chr2H"&cluster==5& PC1 <= -0.005 , cluster:=2]


Repcafor2(pq1,"chr2H",cov,v,2,4,0) -> pq1

pq1[chr=="chr2H"&cluster==5& PC1 >= 0 , cluster:=4]
pq1[chr=="chr2H"&cluster==5& PC1 <= 0 , cluster:=2]


cal_chrlen(pq1,"chr2H",45)


Repcafor2check(pq1,"chr2H",cov,v,1,2,"chr2Hcheck",-0.005)

Repcafor2check(pq1,"chr2H",cov,v,1,3,"chr2Hcheck",0)

Repcafor2check(pq1,"chr2H",cov,v,1,4,"chr2Hcheck",0.01)

Repcafor2check(pq1,"chr2H",cov,v,2,3,"chr2Hcheck",0)

Repcafor2check(pq1,"chr2H",cov,v,2,4,"chr2Hcheck",0)

Repcafor2check(pq1,"chr2H",cov,v,3,4,"chr2Hcheck",0)


####chr3


Repcafor2(pq1,"chr3H",cov,v,3,4,0) -> pq1

pq1[chr=="chr3H"&cluster==5& PC1 >= 0, cluster:=3]
pq1[chr=="chr3H"&cluster==5& PC1 <= 0, cluster:=4]


Repcafor2check(pq1,"chr3H",cov,v,1,2,"chr3Hcheck",0) 

Repcafor2check(pq1,"chr3H",cov,v,1,3,"chr3Hcheck",0)  

Repcafor2check(pq1,"chr3H",cov,v,1,4,"chr3Hcheck",0.01) 

Repcafor2check(pq1,"chr3H",cov,v,2,3,"chr3Hcheck",0)

Repcafor2check(pq1,"chr3H",cov,v,2,4,"chr3Hcheck",0)

Repcafor2check(pq1,"chr3H",cov,v,3,4,"chr3Hcheck",0)



####chr4 


Repcafor2check(pq1,"chr4H",cov,v,1,2,"chr4Hcheck",0) 

Repcafor2check(pq1,"chr4H",cov,v,1,3,"chr4Hcheck",0)  

Repcafor2check(pq1,"chr4H",cov,v,1,4,"chr4Hcheck",0.02) 

Repcafor2check(pq1,"chr4H",cov,v,2,3,"chr4Hcheck",0)

Repcafor2check(pq1,"chr4H",cov,v,2,4,"chr4Hcheck",0)

Repcafor2check(pq1,"chr4H",cov,v,3,4,"chr4Hcheck",0)




####chr5


Repcafor2check(pq1,"chr5H",cov,v,1,2,"chr5Hcheck",0) 

Repcafor2check(pq1,"chr5H",cov,v,1,3,"chr5Hcheck",0)  

Repcafor2check(pq1,"chr5H",cov,v,1,4,"chr5Hcheck",0) 

Repcafor2check(pq1,"chr5H",cov,v,2,3,"chr5Hcheck",0)

Repcafor2check(pq1,"chr5H",cov,v,2,4,"chr5Hcheck",0)

Repcafor2check(pq1,"chr5H",cov,v,3,4,"chr5Hcheck",0)





###chr6
cal_chrlen(pq1,"chr6H",45)

Repcafor2check(pq1,"chr6H",cov,v,1,2,"chr6Hcheck",0) 

Repcafor2check(pq1,"chr6H",cov,v,1,3,"chr6Hcheck",0)  

Repcafor2check(pq1,"chr6H",cov,v,1,4,"chr6Hcheck",0) 

Repcafor2check(pq1,"chr6H",cov,v,2,3,"chr6Hcheck",0)

Repcafor2check(pq1,"chr6H",cov,v,2,4,"chr6Hcheck",0)

Repcafor2check(pq1,"chr6H",cov,v,3,4,"chr6Hcheck",0)



####chr7 
cal_chrlen(pq1,"chr7H",45)


Repcafor2(pq1,"chr7H",cov,v,1,4,0) -> pq1

pq1[chr=="chr7H"&cluster==5& PC1>=0.009, cluster:=1]
pq1[chr=="chr7H"&cluster==5& PC1<=0.009, cluster:=4]



Repcafor2(pq1,"chr7H",cov,v,2,4,0) -> pq1

pq1[chr=="chr7H"&cluster==5& PC1>=0.02, cluster:=2]
pq1[chr=="chr7H"&cluster==5& PC1<=0.02, cluster:=4]


Repcafor2check(pq1,"chr7H",cov,v,1,2,"chr7Hcheck",0) 

Repcafor2check(pq1,"chr7H",cov,v,1,3,"chr7Hcheck",0.01)  

Repcafor2check(pq1,"chr7H",cov,v,1,4,"chr7Hcheck",0.009) 

Repcafor2check(pq1,"chr7H",cov,v,2,3,"chr7Hcheck",0)

Repcafor2check(pq1,"chr7H",cov,v,2,4,"chr7Hcheck",0.02)

Repcafor2check(pq1,"chr7H",cov,v,3,4,"chr7Hcheck",0)


saveRDS(pq1, file="A40_assembly_v2_haplotype_separation_HiC_pq2.Rds")

cal_chrlen(pq1,"chr1H",45)
cal_chrlen(pq1,"chr2H",45)
cal_chrlen(pq1,"chr3H",45)
cal_chrlen(pq1,"chr4H",45)
cal_chrlen(pq1,"chr5H",45)
cal_chrlen(pq1,"chr6H",45)
cal_chrlen(pq1,"chr7H",45)


> cal_chrlen(pq1,"chr1H",45)
[1] 394050661
[1] 352238082
[1] 379430333
[1] 350171965
> cal_chrlen(pq1,"chr2H",45)
[1] 445482660
[1] 451085238
[1] 504001382
[1] 465118090
> cal_chrlen(pq1,"chr3H",45)
[1] 456338782
[1] 474787619
[1] 473146615
[1] 488562875
> cal_chrlen(pq1,"chr4H",45)
[1] 381461783
[1] 417103232
[1] 397240082
[1] 388591793
> cal_chrlen(pq1,"chr5H",45)
[1] 447148528
[1] 412256251
[1] 414050167
[1] 418149355
> cal_chrlen(pq1,"chr6H",45)
[1] 374488285
[1] 345488812
[1] 336125826
[1] 370503325
> cal_chrlen(pq1,"chr7H",45)
[1] 438593928
[1] 474483036
[1] 440080902
[1] 450691554


readRDS("A40_assembly_v2_haplotype_separation_HiC_pq2.Rds")->pq2

pq2[cluster==0, cluster := 22]
pq2[, .(contig,cluster)][pp, on="contig"] -> pq4


readRDS("A40_assembly_v2_haplotype_separation_HiC.Rds")->pp

cov[, .(contig=scaffold, cc)][pp, on="contig"] -> pp

pp[cc<45, col := "red"]
pp[cc>=45&cc<80,col:="blue"]
pp[cc>=80&cc<120,col:="black"]
pp[cc>=120,col:="yellow"]

pdf("A40_haplotype_separation_HiC_pc1_2_precheck.pdf", height=10, width=10)

lapply(morexfai[1:7, chr], function(i){
  pp[chr == i, plot(PC2, PC1, col=col, pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC1",
            xlab="PC2")]
  legend("topleft", pch=19, bty='n', legend=c("1X", "2X","3X","4X"), col=c("red", "blue","black","orange"))
 })

dev.off()


pdf("A40_haplotype_separation_HiC_after_cluster_pc1_2.pdf", height=10, width=10)

lapply(morexfai[1:7, chr], function(i){
  pq4[chr == i, plot(PC2, PC1, col=cluster, pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC1",
            xlab="PC2")]
  legend("topleft", pch=19, bty='n', legend=c("1","2","3","4","22","33","44"), col=c(1,2,3,4,22,33,44))
 })

dev.off()

###To identify which haplotypes the contigs with 2X, 3X, and 4X coverage belong to.

Reassign <- function(pq4, i, n)
{
j=pq4[contig==i]$chr
c(pq4[contig==i]$PC1, pq4[contig==i]$PC2, pq4[contig==i]$PC3, pq4[contig==i]$PC4) -> point1
c(mean(pq4[cluster==1 & chr ==j]$PC1),mean(pq4[cluster==1& chr==j]$PC2),mean(pq4[cluster==1& chr==j]$PC3),mean(pq4[cluster==1& chr==j]$PC4)) -> center1
c(mean(pq4[cluster==2 & chr ==j]$PC1),mean(pq4[cluster==2& chr==j]$PC2),mean(pq4[cluster==2& chr==j]$PC3),mean(pq4[cluster==2& chr==j]$PC4)) -> center2
c(mean(pq4[cluster==3 & chr ==j]$PC1),mean(pq4[cluster==3& chr==j]$PC2),mean(pq4[cluster==3& chr==j]$PC3),mean(pq4[cluster==3& chr==j]$PC4)) -> center3
c(mean(pq4[cluster==4 & chr ==j]$PC1),mean(pq4[cluster==4& chr==j]$PC2),mean(pq4[cluster==4& chr==j]$PC3),mean(pq4[cluster==4& chr==j]$PC4)) -> center4
d1 <- as.matrix(dist(rbind(point1,center1), method = "euclidean"))["point1","center1"]
d2 <- as.matrix(dist(rbind(point1,center2), method = "euclidean"))["point1","center2"]
d3 <- as.matrix(dist(rbind(point1,center3), method = "euclidean"))["point1","center3"]
d4 <- as.matrix(dist(rbind(point1,center4), method = "euclidean"))["point1","center4"]
c(d1,d2,d3,d4)-> dd 
sort(dd) -> dd 
cl=""
if(d1<=dd[n]){
     cl=paste(cl,"1",sep="")
}
if(d2<=dd[n]) {
     cl=paste(cl,"2",sep="")
}
if(d3<=dd[n]) {
     cl=paste(cl,"3",sep="")
}
if(d4<=dd[n]) {
     cl=paste(cl,"4",sep="")
}
cl<-as.numeric(cl)
pq4[contig==i,cluster:=cl]
return(pq4)
}

for (i in pq2[cluster==22]$contig){
    Reassign(pq4,i,2)->pq4
}

for (i in pq2[cluster==33]$contig){
    Reassign(pq4,i,3)->pq4
}

setorder(pq4, chr, pos) -> pq4
pq4<-pq4[, idx := 1:.N]

setnames(pq4,"cluster","hap")

#save haplotype separation results
saveRDS(pq4, file="A40_assembly_v2_haplotype_separation_v0.Rds")


#### draw a picture for different haplotypes

pq4 -> pq5

pq5[cc<45, col := "red"]
pq5[cc>=45&cc<80,col:="blue"]
pq5[cc>=80&cc<120,col:="black"]
pq5[cc>=120,col:="orange"]

pq5[hap==12]$hap=5
pq5[hap==13]$hap=6
pq5[hap==14]$hap=7
pq5[hap==23]$hap=8
pq5[hap==24]$hap=9
pq5[hap==34]$hap=10

pq5[hap>=100]$hap=11


#plot results (PC1 score and coverage)
pdf("A40_assembly_v2_haplotype_separation_HiC.pdf", height=8, width=10)
par(mfrow=c(2, 1))
lapply(c("hap"), function(j){
 lapply(morexfai[1:7, chr], function(i){
  pq5[chr == i, plot(pos/1e6, las=1, bty='l', hap, type='n', main=sub("chr", "", chr[1]), ylab="hap",
        xlab="Hv syntenic position [Mb]", xlim=c(0, morexfai[i, len/1e6, on="chr"]), ylim=c(0,12))]
  pq5[chr == i, lines(lwd=3, c(pos/1e6, (pos + contig_length)/1e6), c(hap, hap), col=hap), by=idx ]
  pq5[chr == i, plot(pos/1e6, las=1, bty='l', cc, type='n', main=sub("chr", "", chr[1]), ylab="coverage",
            xlab="Hv syntenic position [Mb]", xlim=c(0, morexfai[i, len/1e6, on="chr"]))]
  pq5[chr == i, lines(lwd=3, c(pos/1e6, (pos + contig_length)/1e6), c(cc, cc), col=col), by=idx]
 })
})
dev.off()

pdf("A40_haplotype_separation_HiC_after_cluster_final12.pdf", height=10, width=10)

lapply(morexfai[1:7, chr], function(i){
  pq5[chr == i, plot(PC2, PC1, col=hap, pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC1", xlab="PC2")]
  legend("topleft", pch=19, bty='n', legend=c("1","2","3","4","5","6","7","8","9","10"), col=c(1,2,3,4,5,6,7,8,9,10))
 })

dev.off()

pdf("A40_haplotype_separation_HiC_after_cluster_final23.pdf", height=10, width=10)

lapply(morexfai[1:7, chr], function(i){
  pq5[chr == i, plot(PC3, PC2, col=hap, pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC2", xlab="PC3")]
  legend("topleft", pch=19, bty='n', legend=c("1","2","3","4","5","6","7","8","9","10"), col=c(1,2,3,4,5,6,7,8,9,10))
 })

dev.off()

pdf("A40_haplotype_separation_HiC_after_cluster_final34.pdf", height=10, width=10)

lapply(morexfai[1:7, chr], function(i){
  pq5[chr == i, plot(PC3, PC4, col=hap, pch=20, las=1,bty='l', type="p", main=sub("chr", "", chr[1]), ylab="PC3", xlab="PC4")]
  legend("topleft", pch=19, bty='n', legend=c("1","2","3","4","5","6","7","8","9","10"), col=c(1,2,3,4,5,6,7,8,9,10))
 })

dev.off()


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
# Assign contigs that are not place in the guide map to haplotype using Hi-C
##

#read haplotype separation and guide map table
readRDS(file="A40_assembly_v2_haplotype_separation_v0.Rds") -> pq
readRDS(file="A40_assembly_v2_Hv_guide+HiClift.Rds") -> cc

readRDS('A40_assembly_v2.Rds') -> assembly_v2
readRDS('A40_assembly_v2_cov.Rds') -> cov


#find and tabulate links between unplaced and placed con tigs
assembly_v2$fpairs[, .(scaffold=scaffold1, pos=pos1, link=scaffold2)] -> ff
assembly_v2$info[, .(scaffold, length)][ff, on="scaffold"] -> ff
ff[pos <= 2e6 | length - pos <= 2e6] -> ff
ff[scaffold != link] -> ff
ff[, .N, key=.(scaffold, link)] -> fa
fa[N > 1] -> fa
fa[scaffold %in% setdiff(cc$scaffold, pq$contig)] -> fa
fa[link %in% pq$contig] -> fa   #### need improve [hap %in% 1:4]
cc[, .(scaffold, scaffold_chr=chr, scaffold_pos=pos)][fa, on="scaffold"] -> fa
cc[, .(link=scaffold, link_chr=chr, link_pos=pos)][fa, on="link"] -> fa
pq[, .(link=contig, hap)][fa, on="link"] -> fa
fa[scaffold_chr == link_chr] -> fa
fa[, .(n=sum(N)), key=.(scaffold, chr=scaffold_chr, hap)] -> fv

fv1 <- fv[hap>=10]

fv1 <- fv1[,id:= paste(scaffold,hap)]

### check one
fv[scaffold=="contig_corrected_v1_3455"]


for (i in fv1$id){
   if(fv1[id==i]$hap==12){
    fv[scaffold == fv1[id==i]$scaffold & hap==1]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==1]$n + fv1[id==i]$n/2
    fv[scaffold == fv1[id==i]$scaffold & hap==2]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==2]$n + fv1[id==i]$n/2
   }
   else if(fv1[id==i]$hap==13){
    fv[scaffold == fv1[id==i]$scaffold & hap==1]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==1]$n + fv1[id==i]$n/2
    fv[scaffold == fv1[id==i]$scaffold & hap==3]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==3]$n + fv1[id==i]$n/2
   }
   else if(fv1[id==i]$hap==14){
    fv[scaffold == fv1[id==i]$scaffold & hap==1]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==1]$n + fv1[id==i]$n/2
    fv[scaffold == fv1[id==i]$scaffold & hap==4]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==4]$n + fv1[id==i]$n/2
   }
   else if(fv1[id==i]$hap==23){
    fv[scaffold == fv1[id==i]$scaffold & hap==2]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==2]$n + fv1[id==i]$n/2
    fv[scaffold == fv1[id==i]$scaffold & hap==3]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==3]$n + fv1[id==i]$n/2
   }
   else if(fv1[id==i]$hap==24){
    fv[scaffold == fv1[id==i]$scaffold & hap==2]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==2]$n + fv1[id==i]$n/2
    fv[scaffold == fv1[id==i]$scaffold & hap==4]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==4]$n + fv1[id==i]$n/2
   }
   else if(fv1[id==i]$hap==34){
    fv[scaffold == fv1[id==i]$scaffold & hap==3]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==3]$n + fv1[id==i]$n/2
    fv[scaffold == fv1[id==i]$scaffold & hap==4]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==4]$n + fv1[id==i]$n/2
   }
   else if(fv1[id==i]$hap==123){
    fv[scaffold == fv1[id==i]$scaffold & hap==1]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==1]$n + fv1[id==i]$n/3
    fv[scaffold == fv1[id==i]$scaffold & hap==2]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==2]$n + fv1[id==i]$n/3
    fv[scaffold == fv1[id==i]$scaffold & hap==3]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==3]$n + fv1[id==i]$n/3
   }
   else if(fv1[id==i]$hap==124){
    fv[scaffold == fv1[id==i]$scaffold & hap==1]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==1]$n + fv1[id==i]$n/3
    fv[scaffold == fv1[id==i]$scaffold & hap==2]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==2]$n + fv1[id==i]$n/3
    fv[scaffold == fv1[id==i]$scaffold & hap==4]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==4]$n + fv1[id==i]$n/3
   }
   else if(fv1[id==i]$hap==134){
    fv[scaffold == fv1[id==i]$scaffold & hap==1]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==1]$n + fv1[id==i]$n/3
    fv[scaffold == fv1[id==i]$scaffold & hap==3]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==3]$n + fv1[id==i]$n/3
    fv[scaffold == fv1[id==i]$scaffold & hap==4]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==4]$n + fv1[id==i]$n/3
   }
   else if(fv1[id==i]$hap==234){
    fv[scaffold == fv1[id==i]$scaffold & hap==2]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==2]$n + fv1[id==i]$n/3
    fv[scaffold == fv1[id==i]$scaffold & hap==3]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==3]$n + fv1[id==i]$n/3
    fv[scaffold == fv1[id==i]$scaffold & hap==4]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==4]$n + fv1[id==i]$n/3
   }
   else if(fv1[id==i]$hap==1234){
    fv[scaffold == fv1[id==i]$scaffold & hap==1]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==1]$n + fv1[id==i]$n/4
    fv[scaffold == fv1[id==i]$scaffold & hap==2]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==2]$n + fv1[id==i]$n/4
    fv[scaffold == fv1[id==i]$scaffold & hap==3]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==3]$n + fv1[id==i]$n/4
    fv[scaffold == fv1[id==i]$scaffold & hap==4]$n <- fv[scaffold == fv1[id==i]$scaffold & hap==4]$n + fv1[id==i]$n/4
   }
}


fv[hap<10] -> fv 
fv[, p := n/sum(n), by=scaffold]
assembly_v2$info[, .(scaffold, length)][fv, on="scaffold"] -> fv
fv[length >= 0][order(-p)][!duplicated(scaffold)][order(-length)] -> hap_lift0

pdf("p_assigment.pdf", height=10, width=10)
plot(density(fv[hap==1]$p),xlim=c(0,1), col="red")
lines(density(fv[hap==2]$p),xlim=c(0,1), col="blue")
lines(density(fv[hap==3]$p),xlim=c(0,1), col="green")
lines(density(fv[hap==4]$p),xlim=c(0,1), col="black")
lines(density(hap_lift0$p),xlim=c(0,1), col="orange")
abline(v=c(0.55, 0.4, 0.3), col="blue")
dev.off()


cov[, .(scaffold, cc)][hap_lift0, on="scaffold"] -> hap_lift0

pdf("p_assigment_coverage.pdf", height=10, width=10)
plot(density(hap_lift0[!is.na(cc)]$cc))
abline(v=c(55, 90), col="blue")
dev.off()




saveRDS(hap_lift0, file="A40_assembly_v2_hap_lift0.Rds") 

readRDS("A40_assembly_v2_hap_lift0.Rds") -> hap_lift0
#keep contigs unassigned to one haplotype OR contigs assigned to both haplotyped with double coverage
hap_lift0 -> hap_lift

#keep onlt contigs >= 300 kb and sum(links) >=50
hap_lift[length >= 3e5] -> hap_lift
hap_lift[n/p >= 50] -> hap_lift


fv[hap==1] -> fv11
setnames(fv11,"p","hap1")
fv[hap==2] -> fv22
setnames(fv22,"p","hap2")
fv[hap==3] -> fv33
setnames(fv33,"p","hap3")
fv[hap==4] -> fv44
setnames(fv44,"p","hap4")

fv11[, .(scaffold,hap1)][fv22, on="scaffold"] -> fv55
fv33[, .(scaffold,hap3)][fv55, on="scaffold"] -> fv55
fv44[, .(scaffold,hap4)][fv55, on="scaffold"] -> fv55

fv66<-data.frame(fv55$scaffold,fv55$hap1,fv55$hap2,fv55$hap3,fv55$hap4)
fv66[is.na(fv66)]=0
data.table(fv66) -> fv66
colnames(fv66)<-c("scaffold","hap1","hap2","hap3","hap4")
fv66[,ma1:=0]
fv66[,ma2:=0]
fv66[,ma3:=0]
fv66[,ma4:=0]
for (i in row(fv66)){
  dd=sort(fv66[i,2:5])
  colnames(dd)<-c("ma1","ma2","ma3","ma4")
  fv66[i]$ma1=dd$ma1
  fv66[i]$ma2=dd$ma2
  fv66[i]$ma3=dd$ma3
  fv66[i]$ma4=dd$ma4
}

pdf("ppfv_ma.pdf", height=10, width=10)
plot(density(fv66$ma4-fv66$ma3),xlim=c(0,1),col="orange")
plot(density(fv66$ma3-fv66$ma2),xlim=c(0,1),col="blue")
plot(density(fv66$ma2-fv66$ma1),xlim=c(0,1),col="green")
dev.off()

fv66[,hap:=0]
fv66[ma4-ma3>=0.3,hap:=11]
fv66[ma4-ma3<0.2 & ma3-ma2 >=0.2,hap:=22]
fv66[ma4-ma3<0.2 & ma3-ma2 <0.2 & ma2-ma1>=0.2, hap:=33]
fv66[ma4-ma3<0.2 & ma3-ma2 <0.2 & ma2-ma1<0.2 & ma1>=0.2, hap:=44]

### Distinguish the haplotype affiliation of small fragments based on phased large fragments and Hi-C links

Reassign_samll <- function(fv, i, n)
{
d1 <- fv[scaffold == i & hap == 1 ]$p
d2 <- fv[scaffold == i & hap == 2 ]$p
d3 <- fv[scaffold == i & hap == 3 ]$p
d4 <- fv[scaffold == i & hap == 4 ]$p
if (length(fv[scaffold == i & hap == 1 ]$p)==0) {
   d1=0
}
if (length(fv[scaffold == i & hap == 2 ]$p)==0) {
   d2=0
}
if (length(fv[scaffold == i & hap == 3 ]$p)==0) {
   d3=0
}
if (length(fv[scaffold == i & hap == 4 ]$p)==0) {
   d4=0
}
c(d1,d2,d3,d4)-> dd 
sort(dd,decreasing=TRUE) -> dd 
cl=""
if(d1>=dd[n]){
     cl=paste(cl,"1",sep="")
}
if(d2>=dd[n]) {
     cl=paste(cl,"2",sep="")
}
if(d3>=dd[n]) {
     cl=paste(cl,"3",sep="")
}
if(d4>=dd[n]) {
     cl=paste(cl,"4",sep="")
}
cl<-as.numeric(cl)
return(cl)
}


for (i in fv66[hap==22]$scaffold){
    Reassign_samll(fv,i,2)->hap_lift[scaffold==i]$hap
}

for (i in fv66[hap==33]$scaffold){
    Reassign_samll(fv,i,3)->hap_lift[scaffold==i]$hap
}

for (i in fv66[hap==44]$scaffold){
    Reassign_samll(fv,i,4)->hap_lift[scaffold==i]$hap
}

for (i in fv66[hap==0]$scaffold){
    hap_lift[scaffold==i]$hap=0
}

hap_lift[hap>=10&hap<=100&cc<45]$hap <- 0
hap_lift[hap>=100&hap<=1000&cc<80]$hap <- 0
hap_lift[hap>=1000&cc<120]$hap<-0 

sum(hap_lift[hap==0]$length)


hap_lift[hap!=0] -> hap_lift

#combine both haplotype tables
rbind(pq[, .(scaffold=contig, hap)],  hap_lift[, .(scaffold, hap)]) -> hh
hh[cc, on="scaffold"] -> hh
saveRDS(hh, file="A40_assembly_v2_Hv_guide+HiClift+hap1.Rds") 

sum(hh[is.na(hap)]$length)

hh -> hh1 

#cov[, .(scaffold, cc)][hh1, on="scaffold"] -> hh1
fread('/filer-dg/agruppen/seq_shared/mascher/morexV3_pseudomolecules_200421/200416_MorexV3_pseudomolecules.fasta.fai', sel=1:2, col.names=c("chr", "len"))->morexfai

setorder(hh1, chr, pos) -> hh1
hh1<-hh1[, idx := 1:.N]


hh1[hap==12]$hap<-5
hh1[hap==13]$hap<-6
hh1[hap==14]$hap<-7
hh1[hap==23]$hap<-8
hh1[hap==24]$hap<-9
hh1[hap==34]$hap<-10

hh1[hap==123]$hap<-11
hh1[hap==124]$hap<-12
hh1[hap==134]$hap<-13
hh1[hap==234]$hap<-14

hh1[hap==1234]$hap<-15

hh1[is.na(hap)]$hap <- 16

pdf("A40_assembly_v2_haplotype_separation_HiC2111.pdf", height=8, width=10)
par(mfrow=c(1, 1))
lapply(c("cluster"), function(j){
 lapply(morexfai[1:7, chr], function(i){
  hh1[chr == i, plot(pos/1e6, las=1, bty='l', hap, type='n', main=sub("chr", "", chr[1]), ylab="hap",
        xlab="Hv syntenic position [Mb]", xlim=c(0, morexfai[i, len/1e6, on="chr"]), ylim=c(0,16))]
  hh1[chr == i, lines(lwd=3, c(pos/1e6, (pos + length)/1e6), c(hap, hap), col=hap), by=idx ]
 })
})

dev.off()

