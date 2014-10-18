library(ggplot2)
library("grid")
require(mgcv)
require(splines)
require(MASS)

prefix ="ceu.yri"
output.name = paste(prefix,"4",sep=".")
data.dir = "K:\\Data\\IBDAdmixed\\New5\\"
seglen = "3.cM"

data = c()

curr.data = read.table(paste(data.dir,output.name,".ibdadmixed.txt.",seglen,".stats.txt",sep=""))
curr.data$Method = "IBDAdmixed"
data = rbind(data,curr.data)

# curr.data = read.table(paste(data.dir,output.name,".ibdadmixed.0.5.txt.stats.txt",sep=""))
# curr.data$Method = "IBDAdmixed>0.5"
# data = rbind(data,curr.data)

# curr.data = read.table(paste(data.dir,output.name,".naive.ibdadmixed.txt.",seglen,".stats.txt",sep=""))
# curr.data$Method = "Naive"
# data = rbind(data,curr.data)

curr.data = read.table(paste(data.dir,output.name,".genos.beagle3.ibd.txt.",seglen,".stats.txt",sep=""))
curr.data$Method = "Beagle3"
data = rbind(data,curr.data)

curr.data = read.table(paste(data.dir,prefix,".parente.ibd.txt.",seglen,".stats.txt",sep=""))
curr.data$Method = "Parente"
data = rbind(data,curr.data)

# curr.data = read.table(paste(data.dir,output.name,".germline.ibd.txt.2cM.stats.txt",sep=""))
# curr.data$Method = "GERMLINE"
# data = rbind(data,curr.data)

colnames(data) = c("Score","Sensitivity","FDR","FPR","Method")
data = data[data$FPR>0.0001,]
cbPalette <- c("#56B4E9", "#0072B2", "#009E73", "#E69F00", "#D55E00", "#CC79A7")

ggplot() + 
  stat_smooth(data=data[data$Method == "IBDAdmixed",], aes(x=FPR, y=Sensitivity, group=Method, colour=Method, shape=Method),
              method = "lm", formula = y ~ ns(x, 5), n=100, size = 1.3, se=FALSE, fullrange=TRUE) + 
  stat_smooth(data=data[data$Method == "Beagle3",], aes(x=FPR, y=Sensitivity, group=Method, colour=Method, shape=Method),
              method = "lm", formula = y ~ ns(x, 6), n=300, size = 1.3, se=FALSE, fullrange=TRUE) + 
  stat_smooth(data=data[data$Method == "Naive",], aes(x=FPR, y=Sensitivity, group=Method, colour=Method, shape=Method),
              method = "lm", formula = y ~ ns(x, 5), n=300, size = 1.3, se=FALSE, fullrange=TRUE) + 
  stat_smooth(data=data[data$Method == "Parente",], aes(x=FPR, y=Sensitivity, group=Method, colour=Method, shape=Method),
              method = "lm", formula = y ~ ns(x, 4), n=300, size = 1.3, se=FALSE, fullrange=TRUE) + 
#   stat_smooth(data=data[data$Method != "IBDAdmixed",], aes(x=FPR, y=Sensitivity, group=Method, colour=Method, shape=Method),
#               method = "lm", formula = y ~ ns(x, 5), n=300, size = 1.3, se=FALSE, fullrange=TRUE) + 
  #stat_smooth(data=data[data$Method == "IBDAdmixed",], aes(x=FPR, y=Sensitivity, group=Method, colour=Method, shape=Method),
  #            method = "gam", formula = y ~ s(x, k = 7), n=300, size = 1.3, se=FALSE, fullrange=TRUE) + 
  geom_point(data=data, aes(x=FPR, y=Sensitivity, group=Method, colour=Method, shape=Method), size=4) + 
  scale_x_continuous(expand = c(0,0), limits = c(0, 0.03)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  scale_colour_manual(values=cbPalette) + 
  scale_shape_manual(name="Method", values=c(22,21,23,24,25,26)) +
  scale_linetype_discrete(name="Method") +
  xlab("False Positive Rate") + ylab("Sensitivity") +
  theme_bw() + 
  theme(legend.justification=c(1,0), 
        legend.position=c(1,0), 
        legend.text = element_text(size = 10), 
        legend.title = element_blank(), 
        legend.key = element_blank(), 
        legend.background = element_rect(colour = "black",size=0.1), 
        legend.key.height=unit(0.8,"line"), 
        panel.grid.major =  element_line(colour = "grey", size = 0.05, linetype = 'dashed'), 
        panel.grid.minor =  element_line(colour = "grey", size = 0.05, linetype = 'dashed'), 
        panel.border = element_rect(colour = "black",size=1),
        axis.title.x = element_text(size=12),
        axis.text.x  = element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.text.y  = element_text(size=10))