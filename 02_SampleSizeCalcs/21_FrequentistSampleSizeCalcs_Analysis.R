#install.packages("writexl")
library("writexl")
library("ggplot2")
library("fs")
library("tidyr")
library("haven")

getwd()
setwd("C:/Users/tofis/Dropbox (IDinsight)/IE Design for DFEs/Code/03_Output")
outpath = "C:/Users/tofis/Dropbox (IDinsight)/IE Design for DFEs/Code/03_Output"


### Load uninformative results
##############################
twosided = read_dta(path(getwd(),"twosided.dta",ext=""))
twosided[dim(twosided)[2]+1]<-"Two-Sided"

onesided = read_dta(path(getwd(),"onesided.dta",ext=""))
onesided[dim(onesided)[2]+1]<-"One-Sided"

SampPlot <- rbind(twosided,onesided)
names(SampPlot)[15]<- "Type"

summary(SampPlot)

# Power graphs 
Samp_graph <- ggplot(data = SampPlot, aes(y=N, x=alpha,group=Type,color=Type)) + 
  geom_line(size=1) + geom_point(size=1.25) +  
  ylab("Total Sample Size (N)") + xlab("Significance Level (alpha)") +
  scale_y_continuous(limits=c(50, 1200), breaks=seq(100,1200,100)) +
  scale_colour_discrete("Sample Size") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        legend.position = "right", legend.justification = "center",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        legend.key.size = unit(0.3, "cm"),
        legend.key.width = unit(0.5,"cm"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10))

Samp_graph
ggsave(file.path(path(outpath,"TwoVSOnesidedTest.png",ext="")), Samp_graph, height=7, width=7)
