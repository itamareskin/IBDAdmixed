library(data.table)
setwd("~/Data/IBDAdmixed/fromLecs")
beagle = fread("beagle4.ibd.long.txt",header=F)
colnames(beagle) = c("ind1","ind2","start","end","score")
beagle$chr = 
beagle$strand = '+'
true.ibd = fread("AfricanAmericans8.trueibd.long.txt",header=F)
colnames(true.ibd) = c("ind1","ind2","start","end","score")
true.ibd$chr = 1
true.ibd$strand = '+'

beagle$chr = beagle$ind1*1000+beagle$ind2
true.ibd$chr = true.ibd$ind1*1000+true.ibd$ind2

library("IRanges")
library("GenomicRanges")
bed <- with(true.ibd[true.ibd$ind1 == ], GRanges(chr, IRanges(start, end), strand, ind1, ind2, score))
bed2 <- with(beagle, GRanges(chr, IRanges(start, end), strand, ind1, ind2, score))
res = intersect(bed, bed)

subject = RangedData(ranges=IRanges(start=true.ibd$start,end=true.ibd$end), ind1 = true.ibd$ind1, ind2 = true.ibd$ind2)
query = RangedData(ranges=IRanges(start=beagle$start,end=beagle$end), ind1 = beagle$ind1, ind2 = beagle$ind2)


ir = IRanges(start=true.ibd$start,end=true.ibd$end)
ir = split(ir,true.ibd$chr)
cov = coverage(ir)