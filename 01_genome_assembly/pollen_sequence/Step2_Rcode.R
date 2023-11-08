#load code for graphical_genotypes() functions
source('/filer-dg/agruppen/DG/mascher/code_repositories/scripts_tutorials/gbs_pipeline/GBS_R_functions.R')

#read the SNP calling results
 fread('210504_Hb_pollenseq_hap1_samtools_normalized.txt.gz')->z
 setnames(z, c("chr", "pos", "ref", "alt", "qual", "sample", "raw", "dp", "dp_ref", "dp_alt"))
 z[dp_ref + dp_alt >= 0.9 * dp] -> vv
 vv[, gt := as.integer(NA)]
 vv[dp_ref >= 0.9 * (dp_alt + dp_ref), gt := 0]
 vv[dp_alt >= 0.9 * (dp_alt + dp_ref), gt := 2]
 vv[dp_ref >= 0.3 * (dp_alt + dp_ref) & dp_ref <= 0.7 * (dp_alt + dp_ref), gt := 1]
 vv[, snp.id := paste0(chr, ":", pos)]
 vv[sample == 'FB19_011_3_pseudomolecules_v6_hap1.fasta', sample := "Sample_hap1"]
 vv[sample == 'FB19_011_3_pseudomolecules_v6_hap2.fasta', sample := "Sample_hap2"]
 saveRDS(vv, file="210504_Hb_pollenseq_hap1_samtools_normalized.Rds")

#select samples where hap2 is different fron hap1 and the Illumina WGS data of hap1+hap2 (parent) is het
 vv[, .N, key=sample] -> ss
 vv[c('Sample_hap1', 'Sample_hap2'), on="sample"] -> vs
 dcast(vs, snp.id ~ sample, value.var=c("dp", "gt")) -> vs
 vs[dp_Sample_hap2 >= 3 & gt_Sample_hap2 == 2] -> vsf
 vsf[dp_Sample_hap1 >= 3 & gt_Sample_hap1 == 0] -> vsf

#get consensus genotypes in 200 kb bi s
 vv[vsf$snp.id, on="snp.id"] -> vf
 options(scipen=20)
 vf[, bin := pos %/% 2e5 * 2e5]
 vf[!is.na(gt), .N, key=.(chr, bin, gt, sample)] -> vb
 vb[gt %in% c(0,2) & !sample %in% c("Sample_hap1", "Sample_hap2")] -> vb
 vb[, gt := paste0("gt", gt)]
 dcast(vb, chr + bin + sample ~ gt, value.var="N", fill=0) -> b
 b[, n := gt0 + gt2]
 b[, r := gt2 / n]
 b[n >= 20] -> bb #select bins with at least 20 SNPs

####decide chr and sample from which hap

 vf[!is.na(gt), .N, key=.(chr, gt, sample)] -> vb
 vb[gt %in% c(0,2) ] -> vb
 vb[, gt := paste0("gt", gt)]
 dcast(vb, chr + sample ~ gt, value.var="N", fill=0) -> b
 b[, n := gt0 + gt2]
 b[, r := gt2 / n]
 b[r<=0.4,hap:="hap1"]
 b[r>=0.6,hap:="hap2"]

 
###send hap to the every SNP
 b[, .(chr, sample, hap)][vf, on=c("chr","sample")]  -> vf1 
 vf1[hap %in% c("hap1","hap2") & gt %in% c(0,2) & !sample %in% c("Sample_hap1", "Sample_hap2")] -> vf1

#####change hap2 to hap1 __________ -> ++++++++++ / ____+________+ -> ++++_++++++++_

 vf1[hap=="hap1",gtn:=gt]
 vf1[hap=="hap2" & gt==0 ,gtn:=2]
 vf1[hap=="hap2" & gt==2 ,gtn:=0]
 vf1[, gtn := paste0("gtn", gtn)]

#### hap2 now is the recombination or phasing error 

 vf1[!is.na(gt), .N, key=.(chr, gtn, snp.id)] -> vb1
 dcast(vb1, chr + snp.id ~ gtn, value.var="N", fill=0) -> b1
 b1[, n := gtn0 + gtn2 ]
 b1[n>=10,r:=gtn0/n]
 b1[n>=10]-> b2
 #count(b2)
 
 ####based hap information to calulate the error rate
 b2
 ##6028524
 b2[r<=0.2]
 #2220
 b2[r<=0.4]
 #30736
 b2[r<0.5]
 #60934 
#### error rate 0.0003682493 ----> 0.005098429 --> 0.01010762 
#99%


Haplotype 1 
00000000000000
++++++++++++++

Haplotype 2
11111111111111
______________

               Decide hap based
gt             all SNP‘s hap      gtn
++++_+++_____+  6/14  hap1  --->  ++++_+++_____+
__++_+++++++++  3/14  hap1  --->  __++_+++++++++
++++_+++++++++  1/14  hap1  --->  ++++_+++++++++
____+________+ 12/14  hap2  --->  ++++_++++++++_ 
____+_________ 13/14  hap2  --->  ++++_+++++++++
                     Count hap1   44550555444444
                        Total     55555555555555
                                  |   |
                                  |   |
          The percent of hap1     0.8 0
               THE ERROR              *    
                                                 The potential 
                                                  error rate
  The level hap2 is caused by      <  0.5 --->       1/14
recombination and phasing error    <= 0.4 --->       1/14
                                   <= 0.2 --->       1/14



#call genotypes in bins
 bb[r < 0.1 | r > 0.9] -> bb
#nrow(unique(bb[, .(chr, bin)]))
#15329
 bb[, gt := as.integer(NA)]
 bb[r < 0.1, gt := 0]
 bb[r > 0.9, gt := 2]
 bb[, pos := bin]
 bb[, snp.id := paste(sep=":", chr, pos)]

 dcast(bb[, .N, key=.(gt=paste0("gt", gt), sample)], sample ~ gt, value.var="N", fill=0) -> sample_summary
 sample_summary[,all_gt:=gt0+gt2]
 keep_list<-sample_summary[all_gt>=10000]$sample

 dcast(bb[, .N, key=.(gt=paste0("gt", gt), chr, pos, snp.id)], snp.id + chr + pos ~ gt, value.var="N")[, n := gt0 + gt2] -> si
 si[, af := gt0/n]
 si[, maf := af]
 si[maf > 0.5, maf := 1 - maf]
 si[maf > 0.1] -> si_good #select bins where MAF > 0.1 (i.e. no extreme segregation distortion)

 nrow(si)
 nrow(si_good)

>  nrow(si)
[1] 15329
>  nrow(si_good)
[1] 14910


 #read reference chromosome lengths
 fread('/filer-dg/agruppen/dg2/fengj/hbulbosum_genome/hb_gbs_analysis/GBSdata/ldhat/pollenseq/FB19_011_3_pseudomolecules_v6_hap1.fasta.fai', sel=1:2)->off
 setnames(off, c("chr", "length"))
 off[, truechr := T]
 rbind(off, data.table(chr="chrUn", truechr=F, length=0)) -> off
 off[, plot_offset := cumsum(c(0, length[-.N] + 1e8))]

 #colors=data.table(gt=0:2, color=c("#0000FF11", "#00FF0011", "#FF000011"))
 #colors=data.table(gt=0:2, color=c("#0000FF33", "#00FF0033", "#FF000033"))
 colors=data.table(gt=0:2, color=c("#FF9326","#A945FF","#0089B2"))
  
 #plot graphical genotypes and filter low coverage sample
 bb[si_good$snp.id, on="snp.id"][sample%in%keep_list] -> y
 y[, sample := sub("Sample_", "", sample)]
 y[, .N, key=sample][order(N)][, idx := 1:.N][, .(sample, idx)] -> sample_order
 graphical_genotypes(t=y, chrlen=off, dist=1e8, wide=F, colors=colors, sample_order=sample_order,
		      margins=c(5,4,4,1), samples_per_page=46, label_dist=400,
 file="210606_Hb_FB19_011_pollen_graphical_genotypes.pdf", title="Graphical genotypes of FB19_011 pollen nuclei")
  #####New
 graphical_genotypes(t=y, chrlen=off, dist=1e8, wide=F, colors=colors, sample_order=sample_order,height=15,
              margins=c(5,4,4,1), samples_per_page=80, label_dist=400,
 file="FB19_011_3_filter_pollen_graphical_genotypes.pdf", title="Graphical genotypes of FB19-011-3 pollen nuclei")

 #plot allele frequency of hap1 along the genome to check for segregation distortion and missing data
 copy(bb)-> good
 good[, .N, key=sample][order(N)][, idx := 1:.N][N >= 10000] -> ss
 good[ss$sample, on="sample"] -> good
 dcast(good[, .N, key=.(gt=paste0("gt", gt), chr, pos, snp.id)], snp.id + chr + pos ~ gt, value.var="N")[, n := gt0 + gt2] -> si
 si[, af := gt0/n]
 si[, maf := af]
 si[maf > 0.5, maf := 1 - maf]

 pdf("210505_hap1_AF.pdf", width=8, height=4)
 par(mfrow=c(2,4))
 lapply(sort(unique(off[chr!="chrUn"]$chr)), function(i){
  si[i, on="chr"][, plot(pch=20, pos/1e6, af*100, las=1, col="#00000033", bty='l', xlab="position in pseudomolecule [Mb]", 
			 ylab="hap1 AF [%]", main=sub("chr", "", i),
	ylim=c(0, 100),
	xlim=c(0, off[i, on="chr"]$length/1e6))]
  abline(col='red', lty=2, h=50)
  abline(col='blue', lty=2, h=85)
  abline(col='blue', lty=2, h=15)
 })
 dev.off()

 pdf("210505_hap1_missing_rates.pdf", width=8, height=4)
 par(mfrow=c(2,4))
 lapply(sort(unique(off[chr!="chrUn"]$chr)), function(i){
  si[i, on="chr"][, plot(pch=20, pos/1e6, n/nrow(ss) * 100, las=1, col="#00000011", bty='l', xlab="position in pseudomolecule [Mb]", 
			 ylab="missing rate [%]", main=sub("chr", "", i),
	ylim=c(0, 100),
	xlim=c(0, off[i, on="chr"]$length/1e6))]
  abline(col='red', lty=2, h=85)
 })
 dev.off()

 #nrow(ss)
 80

 #nrow(si)
 #15694
 si[maf > 0.15 & n >= nrow(ss) * 0.9] -> si_good
 #nrow(si_good)
 #11983
 #si[n < nrow(ss) * 0.9]

 #make genetic map with ASMAP/Mstmap
 good[si_good$snp.id, on="snp.id"][, gt2 := NULL] -> good2
 makeMap(g=good2, pop.type="DH", assign_lg=F, biparental=F, p.value=1e-6)->good_map

 good_map[, nlg := .N, by=lg]
 good_map[good_map[nlg >= 50][, .N, key=.(lg, chr)][, .(lg, chr)], on=c("chr", "lg")] -> zz
 zz[, .(ori=cor(cM, pos, method='s'), max_cM=max(cM)), key=lg][zz, on="lg"] -> zz
 zz[ori > 0, ppos := cM]
 zz[ori < 0, ppos := max_cM - cM]

 #manually define boundaries of recombing regions
 rbind(
 data.table(chr="chr1H", arm="S", pos=50),
 data.table(chr="chr1H", arm="L", pos=410),
 data.table(chr="chr2H", arm="S", pos=50),
 data.table(chr="chr2H", arm="L", pos=520),
 data.table(chr="chr3H", arm="S", pos=35),
 data.table(chr="chr3H", arm="L", pos=480),
 data.table(chr="chr4H", arm="S", pos=30),
 data.table(chr="chr4H", arm="L", pos=430),
 data.table(chr="chr5H", arm="S", pos=30),
 data.table(chr="chr5H", arm="L", pos=380),
 data.table(chr="chr6H", arm="S", pos=50),
 data.table(chr="chr6H", arm="L", pos=410),
 data.table(chr="chr7H", arm="S", pos=60),
 data.table(chr="chr7H", arm="L", pos=490)
 ) -> rb
 saveRDS(rb, file="210507_Hb_FB19_011_map_rec_boundaries.Rds")

 #read barley boundaries
 readRDS('/filer-dg/agruppen/seq_shared/mascher/hordeum_bulbosum_fb19_011_ccs_assembly_201111/Hb_FB19_011_hifiasm_201217/210607_Hb_Hv_POPSEQ_boundaries.Rds')->pb



zz[,.(chr,pos,cM)]->recom
recom[,pos1:=round(pos/10000000,0)]

recom[,.(max(cM)-Min(cM)),key=.(chr,pos1)]
recom[,.(max(cM)-min(cM)),key=.(chr,pos1)]-> recom2

 pdf("RE_no_bar.pdf", height=4, width=8)
 par(mfrow=c(2,4))
 lapply(sort(unique(off[chr!="chrUn"]$chr)), function(i){
  recom2[i, on="chr"][, plot(pch=20, pos1*10, V1/10, las=1, bty='l', xlab="position in pseudomolecule [Mb]", 
             ylab="recombination [cM]", main=sub("chr", "", i),
    ylim=c(0, max(recom2$V1)),
    xlim=c(0, off[i, on="chr"]$length/1e6))]
 # rb[chr == i][, abline(lty=2, lwd=2, col=2, v=pos)]
 # pb[chr == i & hap == 1][, abline(lty=2, lwd=2, col="blue", v=p/1e6)]
 })
 dev.off()

 #plot genetic against physical map
 pdf("210505_LG_chr_correlation_no_bar.pdf", height=4, width=8)
 par(mfrow=c(2,4))
 lapply(sort(unique(off[chr!="chrUn"]$chr)), function(i){
  zz[i, on="chr"][, plot(pch=20, pos/1e6, ppos, las=1, bty='l', xlab="position in pseudomolecule [Mb]", 
			 ylab="position in genetic map [cM]", main=sub("chr", "", i),
	ylim=c(0, unique(zz[, .(max_cM, chr)])[i, on="chr"]$max_cM),
	xlim=c(0, off[i, on="chr"]$length/1e6))]
 # rb[chr == i][, abline(lty=2, lwd=2, col=2, v=pos)]
 # pb[chr == i & hap == 1][, abline(lty=2, lwd=2, col="blue", v=p/1e6)]
 })
 dev.off()

 #plot again with boundaries
 pdf("210505_LG_chr_correlation.pdf", height=4, width=8)
 par(mfrow=c(2,4))
 lapply(sort(unique(off[chr!="chrUn"]$chr)), function(i){
  zz[i, on="chr"][, plot(pch=20, pos/1e6, ppos, las=1, bty='l', xlab="position in pseudomolecule [Mb]", 
			 ylab="position in genetic map [cM]", main=sub("chr", "", i),
	ylim=c(0, unique(zz[, .(max_cM, chr)])[i, on="chr"]$max_cM),
	xlim=c(0, off[i, on="chr"]$length/1e6))]
  rb[chr == i][, abline(lty=2, lwd=2, col=2, v=pos)]
  pb[chr == i & hap == 1][, abline(lty=2, lwd=2, col="blue", v=p/1e6)]
 })
 dev.off()

 #one plot per page
 pdf("210505_LG_chr_correlation_single.pdf", height=4, width=12)
 lapply(sort(unique(off[chr!="chrUn"]$chr)), function(i){
  zz[i, on="chr"][, plot(pch=20, pos/1e6, ppos, las=1, bty='l', xlab="position in pseudomolecule [Mb]", 
			 ylab="position in genetic map [cM]", main=sub("chr", "", i),
	ylim=c(0, unique(zz[, .(max_cM, chr)])[i, on="chr"]$max_cM),
	xlim=c(0, off[i, on="chr"]$length/1e6))]
  rb[chr == i][, abline(lty=2, lwd=2, col=2, v=pos)]
 })
 dev.off()

 #plot for 4H
 fread('/filer-dg/agruppen/seq_shared/mascher/hordeum_bulbosum_fb19_011_ccs_assembly_201111/Hb_FB19_011_hifiasm_201217/./210211_hic_map_v1_hap1_pseudomolecule_AGP.bed')->a
 a[V1 == "chr4H"] -> a

 pdf("210505_LG_chr_correlation_4H.pdf", height=4, width=4)
 i="chr4H"
  zz[i, on="chr"][, plot(pch=20, pos/1e6, ppos, las=1, bty='l', xlab="position in pseudomolecule [Mb]", 
			 ylab="position in genetic map [cM]", main=sub("chr", "", i),
	ylim=c(0, unique(zz[, .(max_cM, chr)])[i, on="chr"]$max_cM),
	xlim=c(400, off[i, on="chr"]$length/1e6))]
 a[, abline(col="gray", v=a$V2/1e6)]
 dev.off()
