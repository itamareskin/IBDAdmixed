data = read.table("~/Dropbox/Computational Biology/IBDAdmixed/figures/AfricanAmericans.artifical.roc.txt")
colnames(data) = c("Method","Score","Sensitivity","FPR")
data = data[data$Method %in% c("IBDAdmixed","fastIBD","PARENTE","IBDAdmixed-Single"),]
library(ggplot2)
library("grid")
#data = data[data$Method == "fastIBD",]
cbPalette <- c("#56B4E9", "#0072B2", "#009E73", "#E69F00", "#D55E00", "#CC79A7")
ggplot(data=data, aes(x=FPR, y=Sensitivity, group=Method, colour=Method, shape=Method)) + 
  geom_line(size=1) + 
  geom_point(size=2.5, fill="white") + 
  xlim(0, 0.002) + 
  ylim(0, 1) +
  scale_x_continuous(expand = c(0,0), limits = c(0,0.001)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  scale_colour_hue(name="Method", l=30)  +
  scale_colour_manual(values=cbPalette) + 
  scale_shape_manual(name="Method", values=c(22,21,23,24)) +
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
        panel.grid.major =  theme_line(colour = "grey", size = 0.05, linetype = 'dashed'), 
        panel.grid.minor =  theme_line(colour = "grey", size = 0.05, linetype = 'dashed'), 
        panel.border = element_rect(colour = "black",size=1),
        #axis.line = theme_segment(colour = "black"),
        #panel.border = theme_blank(),
        axis.title.x = element_text(size=12),
        axis.text.x  = element_text(size=10),
        axis.title.y = element_text(size=12),
        axis.text.y  = element_text(size=10))

# b <- ggplot(data=data, aes(x=FPR, y=Sensitivity, group=Method, colour=Method, shape=Method)) + 
#   geom_line(size=1) + 
#   geom_point(size=2.5, fill="white") + 
#   xlim(0, 0.002) +    
#   ylim(0, 1) +
#   scale_x_continuous(expand = c(0,0), limits = c(0,0.0025)) +
#   scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +  
#   scale_colour_hue(name="Method", l=30)  +     
#   scale_colour_manual(values=cbPalette) + 
#   scale_shape_manual(name="Method", values=c(22,21,23,24)) +
#   scale_linetype_discrete(name="Method") +
#   xlab("") + ylab("") +
#   theme_bw() + 
#   theme(legend.position="none",
#         panel.grid.major =  theme_line(colour = "grey", size = 0.05, linetype = 'dashed'),
#         panel.grid.minor =  theme_line(colour = "grey", size = 0.05, linetype = 'dashed'),
#         #axis.line = theme_segment(colour = "black"),
#         #panel.border = theme_blank(),
#         axis.title.x = element_text(size=12),
#         axis.text.x  = element_text(size=10),
#         axis.title.y = element_text(size=12),
#         axis.text.y  = element_text(size=10))

# library(gridExtra)
# vpb_ <- viewport(width = 1, height = 1, x = 0.5, y = 0.5)  # the larger map
# vpa_ <- viewport(width = 0.4, height = 0.4, x = 0.7, y = 0.35)  # the inset in upper right
# print(a, vp = vpb_)
#print(b, vp = vpa_)