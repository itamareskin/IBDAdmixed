
library(data.table)
library(reshape)
dir = '/home/eskin/Data/IBDAdmixed/fromLecs'
ibd.admixed = fread(paste(dir,'AfricanAmericans8.IBDAdmixed3.windows.dat',sep="/"))
setnames(ibd.admixed, colnames(ibd.admixed), c("ind1","ind2",seq(length(colnames(ibd.admixed))-2)))
ibd.admixed <- melt(ibd.admixed, id=c("ind1","ind2"))
setnames(ibd.admixed, colnames(ibd.admixed), c("ind1","ind2","window","ibd"))

beagle = fread(paste(dir,'beagle.ibd.windows.filt.txt',sep='/'))
setnames(beagle, colnames(beagle), c("ind1","ind2",seq(length(colnames(beagle))-2)))
beagle = beagle[]
beagle <- melt(beagle, id=c("ind1","ind2"))
setnames(beagle, colnames(beagle), c("ind1","ind2","window","ibd"))

beagle = fread(paste(dir,'beagle.ibd.windows.long.filt.txt',sep='/'))
setnames(beagle, colnames(beagle), c("ind1","ind2","window","ibd","score"))

true.ibd = fread(paste(dir,'AfricanAmericans8.trueibd.windows.filt.dat',sep='/'))
setnames(true.ibd, colnames(true.ibd), c("ind1","ind2",seq(length(colnames(true.ibd))-2)))
true.ibd <- melt(true.ibd, id=c("ind1","ind2"))
setnames(true.ibd, colnames(true.ibd), c("ind1","ind2","window","ibd"))

true.ibd.len = fread(paste(dir,'AfricanAmericans8.trueibd.windows.filt.dat.len.dat',sep='/'))
setnames(true.ibd.len, colnames(true.ibd.len), c("ind1","ind2",seq(length(colnames(true.ibd.len))-2)))
true.ibd.len <- melt(true.ibd.len, id=c("ind1","ind2"))
setnames(true.ibd.len, colnames(true.ibd.len), c("ind1","ind2","window","ibd.len"))


data = merge(ibd.admixed,beagle,by=c("ind1","ind2","window"),all.x=T)
data = merge(data,true.ibd,by=c("ind1","ind2","window"),all.x=T)
data = merge(data,true.ibd.len,by=c("ind1","ind2","window"),all.x=T)
setnames(data, colnames(data), c("ind1","ind2","window","ibd.admixed","beagle","beagle.score","true.ibd","true.ibd.len"))
data$beagle[is.na(data$beagle)] = 0
data$true.ibd[is.na(data$true.ibd)] = 0

data.new = data

windows = dim(data.new)[1]
beagle.TP = sum(data.new$true.ibd == 1 & data.new$beagle == 1)
ibd.admixed.TP = sum(data.new$true.ibd == 1 & data.new$ibd.admixed == 1)
beagle.FP = sum(data.new$true.ibd == 0 & data.new$beagle == 1)
ibd.admixed.FP = sum(data.new$true.ibd == 0 & data.new$ibd.admixed == 1)
beagle.TN = sum(data.new$true.ibd == 0 & data.new$beagle == 0)
ibd.admixed.TN = sum(data.new$true.ibd == 0 & data.new$ibd.admixed == 0)
beagle.FN = sum(data.new$true.ibd == 1 & data.new$beagle == 0)
ibd.admixed.FN = sum(data.new$true.ibd == 1 & data.new$ibd.admixed == 0)
beagle.sensitivity = beagle.TP/(beagle.TP+beagle.FN)
ibd.admixed.sensitivity = ibd.admixed.TP/(ibd.admixed.TP+ibd.admixed.FN)
beagle.specificity = beagle.TN/(beagle.TN+beagle.FP)
ibd.admixed.specificity = ibd.admixed.TN/(ibd.admixed.TN+ibd.admixed.FP)

roc = data.table(matrix(c(beagle.sensitivity,1-beagle.specificity,ibd.admixed.sensitivity,1-ibd.admixed.specificity),nrow=2,ncol=2,byrow=T))
colnames(roc) = c("Power","FP")
roc$method = c("beagle","ibd.admixed")
library(ggplot2)
p <- ggplot(roc, aes(FP, Power, shape=method))
p + geom_point(size=5) + xlim(0,1) + ylim(0,1) + theme_bw() + opts(legend.text = element_text(size = 16, face = "bold"))


data = read.table("~/Data/IBDAdmixed/fromLecs/artificial3.parente.txt")
colnames(data) = c("score","sensitivity (\\%)","FPR (\\%)")
data$sensitivity = data$sensitivity*100
data$FPR = data$FPR*100
plot(data$FPR, data$sensitivity)
xtable(data,digits=5,)

data = read.table("~/Dropbox/Computational Biology/IBDAdmixed/figures/AfricanAmericans.regularFPR.txt")
colnames(data) = c("method","sensitivity","FDR","FPR")
data$sensitivity = data$sensitivity*100
data$FDR = data$FDR*100
data$FPR = data$FPR*100
colnames(data) = c("method","sensitivity (\\%)","FDR (\\%)","FPR (\\%)")
xtable(data,digits=5,caption="African Americans - default parameters (FPR not controlled)")