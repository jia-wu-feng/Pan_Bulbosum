#!/bin/zsh

tabix='/opt/Bio/bcftools/1.10.2/bin/tabix'
bgzip='/opt/Bio/bcftools/1.10.2/bin/bgzip'
bcftools='/opt/Bio/bcftools/1.10.2/bin/bcftools'

awkcmd='BEGIN{
 OFS=FS="\t"
 if(header){
  print "##fileformat=VCFv4.0"
  print "##FILTER=<ID=PASS,Description=\"All filters passed\">"
  print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
  print "##FORMAT=<ID=DP,Number=.,Type=Integer,Description=\"Read depth\">"
  print "##FORMAT=<ID=DV,Number=.,Type=Integer,Description=\"Read depth of the alternative allele\">"
 } else {
  h["0/0"]=0
  h["0/1"]=1
  h["1/1"]=2
  if(!dphom)
   dphom=1
  if(!dphet)
   dphet=1
  if(!tol)
   tol=0.2499
  if(bed){
   while(getline < bed){
    chr[$1]=$4
    offset[$1]=$5
   }
  }
 }
}

header && /^#(#reference=|CHROM)/

header && /^[^#]/ {
 exit
}

/^#/ || $4 == "N" {
 next
}

{
 o=""
 n=A=B=H=0

 for(i = 10; i <= NF; i++){
  split($i, a, ":")
  split(a[3], b, ",")
  dp = 0+b[1]+b[2]
  dv = 0+b[2]
  if(dp > 0)
   r=dv/dp
  else if($5 == "."){
   o = o"\t0/0:"dp":"dv
   A++
  }
  if(dp >= dphom && r <= tol ){
   o = o"\t0/0:"dp":"dv
   A++
  }
  else if(dp >= dphom && r >= 1-tol){
   o = o"\t1/1:"dp":"dv
   B++
  }
  else if(dp >= dphet && r >= 0.5-tol && r <= 0.5+tol){
   o = o"\t0/1:"dp":"dv
   H++
  }
  else{
   o = o"\t./.:"dp":"dv
   n++
  }
 }
 present = A + B + H
}

!present || A < minhomn || B < minhomn || present < minpresent * (present + n) || A + B < minhomp * present {
 next
}

{
 if(B > A)
  m = (2*A + H) / 2 / present
 else
  m = (2*B + H) / 2 / present
}

m >= minmaf {
 if(bed){
  print chr[$1]"\t"offset[$1]+$2"\t"
 } else {
  printf $1"\t"$2"\t"
 }
 print $1":"$2":"$4":"$5, $4, $5, $6, $7, ".", "GT:DP:DV"o
}'

echo $HOST 1>&2
echo $SLURM_JOB_ID 1>&2
echo $* 1>&2
date 1>&2

zparseopts -D -K -- -out:=o -vcf:=v -dphom:=d -dphet:=e \
 -minqual:=q -mindp:=md -minhomp:=mino -minhomn:=mine -tol:=t \
 -minmaf:=f -minpresent:=p -threads:=tr -bed:=b

bed=$bed[2]
dphom=$d[2]
dphet=$e[2]
minmaf=$f[2]
mindp=$md[2]
minhomn=$mine[2]
minhomp=$mino[2]
out=$o[2]
minpresent=$p[2]
minqual=$q[2]
threads=$tr[2]
tol=$t[2]
list=$v[2]

if [[ -z $bed ]]; then
 bed=0
fi

mktemp | read tmp 
echo "$awkcmd" > $tmp

{
 head -n 1 $list | xargs 2>/dev/null pigz -cd | awk -f $tmp -v header=1 | $bgzip
 parallel -a $list -k -j $threads --will-cite \
  $bcftools view --exclude-types snps -i '"QUAL >='$minqual'"' '{}' \
   \| $bcftools norm -m -any \| $bcftools norm -d indels \| awk -f $tmp \
    -v header=0 \
    -v bed=$bed \
    -v dphom=$dphom \
    -v dphet=$dphet \
    -v minqual=$minqual \
    -v minhomn=$minhomn \
    -v minhomp=$minhomp \
    -v tol=$tol \
    -v minmaf=$minmaf \
    -v minpresent=$minpresent \| $bgzip 
} > $out.vcf.gz && $tabix -Cp vcf $out.vcf.gz

rm -f $tmp

date 1>&2
