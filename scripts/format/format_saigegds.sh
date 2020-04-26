#!/usr/bin/env bash

gunzip -c $1 | awk '
BEGIN {FS=OFS="\t"}
NR==1 {for(i=1;i<=NF;i++) a[$i]=i; print "CHR","POS","Allele1","Allele2","AF_Allele2","imputationInfo","BETA","SE","p.value","N"}
NR >1 {print $a["chr"],$a["pos"],$a["ref"],$a["alt"],$a["AF.alt"],2,$a["beta"],$a["SE"],$a["pval"],$a["num"]}' | \
bgzip -@4 > $1.formatted.gz
