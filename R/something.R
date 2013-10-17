roc = read.table("/home/eskin/Data/IBDAdmixed/fromLecs/parente.roc.txt")
colnames(roc) = c("threshold","sensitivity","FPR")
plot(roc$FPR,roc$sensitivity)

ceu = read.table("/home/eskin/Data/IBDAdmixed/fromLecs/HapMap3_CEU_chr2.map")
colnames(ceu) = c("chrom","rsid","dist","pos")
yri = read.table("/home/eskin/Data/IBDAdmixed/fromLecs/HapMap3_YRI_chr2.map")
colnames(yri) = c("chrom","rsid","dist","pos")

diff = ceu$dist - yri$dist
summary(diff)