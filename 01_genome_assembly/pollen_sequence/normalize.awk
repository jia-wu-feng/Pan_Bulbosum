#!/bin/awk -f

BEGIN{
 OFS="\t"
}

/^#CHROM/ && !header {
 for(i = 1; i <= NF; i++)
  s[i]=$i
 header=1
}

/^#/ { 
 next
}


{
 for(i = 10; i <= NF; i++){
  split($i, a, ":")
  split(a[4], b, ",")
  print $1, $2, $4, $5, $6, s[i], a[1], a[3], b[1], b[2]
 }
}
