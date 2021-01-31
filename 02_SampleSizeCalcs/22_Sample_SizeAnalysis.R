
#install.packages("writexl")
library("writexl")
library("ggplot2")
library("fs")
library("tidyr")
library("ggpubr")

getwd()
setwd("C:/Users/tofis/Dropbox (IDinsight)/IE Design for DFEs/Code/03_Output/01_SampleSizeCalcs")
outpath = "C:/Users/tofis/Dropbox (IDinsight)/IE Design for DFEs/Code/03_Output/01_SampleSizeCalcs"

### Load uninformative results
##############################
load(path(getwd(),"00_Rdata","BayesSampleSizeCalcs_freqCritOpsstat_N1000.RData",ext=""))
write_xlsx(type1err,path(getwd(),"02_Type1Err","type1err.xlsx"))
write_xlsx(freq_power,path(getwd(),"01_Power","power.xlsx"))

# 5000
load(path(getwd(),"00_Rdata","BayesSampleSizeCalcs_freqCritOpsstat_N5000only.RData",ext=""))
write_xlsx(type1err,path(getwd(),"02_Type1Err","type1err_N5000only.xlsx"))
write_xlsx(freq_power,path(getwd(),"01_Power","power_N5000only.xlsx"))

# Load optimist results 
load(path(getwd(),"00_Rdata","BayesSampleSizeCalcs_freqCritOpsstat_N1000_optimisitcprior.RData",ext=""))
write_xlsx(type1err,path(getwd(),"02_Type1Err","type1err_optimistic.xlsx"))
write_xlsx(freq_power,path(getwd(),"01_Power","power_optimistic.xlsx"))

# Load very optimist results 
load(path(getwd(),"00_Rdata","BayesSampleSizeCalcs_freqCritOpsstat_N1000_veryoptimisitcprior.RData",ext=""))
write_xlsx(type1err,path(getwd(),"02_Type1Err","type1err_veryoptimistic.xlsx"))
write_xlsx(freq_power,path(getwd(),"01_Power","power_veryoptimistic.xlsx"))

# Load pessimistic results 
load(path(getwd(),"00_Rdata","BayesSampleSizeCalcs_freqCritOpsstat_N1000_pessimmisticprior.RData",ext=""))
write_xlsx(type1err,path(getwd(),"02_Type1Err","type1err_pessimistic.xlsx"))
write_xlsx(freq_power,path(getwd(),"01_Power","power_pessimistic.xlsx"))

# Load very pessimistic results 
load(path(getwd(),"00_Rdata","BayesSampleSizeCalcs_freqCritOpsstat_N1000_verypessimmisticprior.RData",ext=""))
write_xlsx(type1err,path(getwd(),"02_Type1Err","type1err_verypessimistic.xlsx"))
write_xlsx(freq_power,path(getwd(),"01_Power","power_verypessimistic.xlsx"))


### combine results 
###################


type1err_neutral <- readxl::read_xlsx(path(getwd(),"02_Type1Err","type1err.xlsx",ext=""))
type1err_neutral[,9] <- "Uninformative"
type1err_neutral_5k <- readxl::read_xlsx(path(getwd(),"02_Type1Err","type1err_N5000only.xlsx",ext=""))
type1err_neutral_5k[,9] <- "Uninformative"
type1err_opt <- readxl::read_xlsx(path(getwd(),"02_Type1Err","type1err_optimistic.xlsx",ext=""))
type1err_opt[,9] <- "Optimistic"
type1err_veropt <- readxl::read_xlsx(path(getwd(),"02_Type1Err","type1err_veryoptimistic.xlsx",ext=""))
type1err_veropt[,9] <- "Very Optimistic"
type1err_pes <- readxl::read_xlsx(path(getwd(),"02_Type1Err","type1err_pessimistic.xlsx",ext=""))
type1err_pes[,9] <- "Pessimistic"
type1err_verpes <- readxl::read_xlsx(path(getwd(),"02_Type1Err","type1err_verypessimistic.xlsx",ext=""))
type1err_verpes[,9] <- "Very Pessimistic"

type1err <- rbind(type1err_neutral,type1err_neutral_5k,type1err_opt,type1err_veropt,type1err_pes,type1err_verpes)
names(type1err)[9]<- "Prior"
type1err_bayes <- type1err[,c(1:5,9)]
h11 = rep("Sample Size ",dim(type1err_bayes[,1])[1]) 
h22 = type1err_bayes %>% pull(samp_size)
type1err_bayes[,1] <- tibble(paste0(h11, h22))


colnames(type1err_bayes)=c("SampleSize","50","60","70","80","Prior")
type1err_bayes_long <- pivot_longer(type1err_bayes,cols=2:5,names_to="Threshold",values_to= "Type1err")
type1err_bayes_long = as.data.frame(type1err_bayes_long)
type1err_bayes_long[,3] <- as.numeric(type1err_bayes_long[,3])

summary(type1err_bayes_long)

#power
power_neutral <- readxl::read_xlsx(path(getwd(),"01_Power","power.xlsx",ext=""))
power_neutral[,9] <- "Uninformative"
power_neutral_5k <- readxl::read_xlsx(path(getwd(),"01_Power","power_N5000only.xlsx",et=""))
power_neutral_5k[,9] <- "Uninformative"
power_opt <- readxl::read_xlsx(path(getwd(),"01_Power","power_optimistic.xlsx",ext=""))
power_opt[,9] <- "Optimistic"
power_pes <- readxl::read_xlsx(path(getwd(),"01_Power","power_pessimistic.xlsx",ext=""))
power_pes[,9] <- "Pessimistic"
power_veropt <- readxl::read_xlsx(path(getwd(),"01_Power","power_veryoptimistic.xlsx",ext=""))
power_veropt[,9] <- "Very Optimistic"
power_verpes <- readxl::read_xlsx(path(getwd(),"01_Power","power_verypessimistic.xlsx",ext=""))
power_verpes[,9] <- "Very Pessimistic"

power <- rbind(power_neutral,power_neutral_5k,power_opt,power_pes,power_veropt,power_verpes)
names(power)[9]<- "Prior"
power_bayes <- power[,c(1:5,9)]
h1 = rep("Sample Size ",dim(power_bayes[,1])[1]) 
h2 = power_bayes %>% pull(samp_size)
power_bayes[,1] <- tibble(paste0(h1, h2))

colnames(power_bayes)=c("SampleSize","50","60","70","80","Prior")
power_bayes_long <- pivot_longer(power_bayes,cols=2:5,names_to="Threshold",values_to= "Power")
power_bayes_long = as.data.frame(power_bayes_long)
power_bayes_long[,3] <- as.numeric(power_bayes_long[,3])

summary(power_bayes_long)

# Power graphs 
Power_graph <- ggplot(data = transform(power_bayes_long,SampleSize=factor(SampleSize,levels=c("Sample Size 50","Sample Size 100","Sample Size 200","Sample Size 300","Sample Size 500","Sample Size 1000","Sample Size 2000","Sample Size 4000","Sample Size 5000"))), aes(x=Threshold, y=Power)) +
  geom_point( aes(shape=Prior,color=Prior), size=0.75, alpha=.8) +
  ylim(0, 1) +  coord_flip() + xlab("Sufficiency Threshold") +
  scale_y_continuous(limits=c(0, 1), breaks=seq(0,1,0.1)) +
  scale_shape_manual(values=c(1,2,15,16,17))+
  scale_color_manual(values=c('#999999','#E69F00', '#56B4E9', '#D9717D','#BECA55')) +
  geom_hline(yintercept=c(0.8, 0.9),color="gray",alpha=0.7) + 
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        legend.position = "right", legend.justification = "center",
        legend.title = element_text(size=5),
        legend.text = element_text(size=5),
        legend.key.size = unit(0.25, "cm"),
        legend.key.width = unit(0.5,"cm"),
        axis.text=element_text(size=4),
        axis.title=element_text(size=5),
        strip.text.x = element_text(size = 4),
        strip.text.y = element_text(size = 4)) +
  facet_wrap(~SampleSize,ncol=1)

Power_graph

# Type 1 err graphs 
summary(type1err_bayes_long)

Type1err_graph <- ggplot(data = transform(type1err_bayes_long,SampleSize=factor(SampleSize,levels=c("Sample Size 50","Sample Size 100","Sample Size 200","Sample Size 300","Sample Size 500","Sample Size 1000","Sample Size 2000","Sample Size 4000","Sample Size 5000"))), aes(x=Threshold, y=Type1err)) + 
  geom_point( aes(shape=Prior,color=Prior), size=0.75, alpha=.8) +
  ylim(0, 0.5) +  coord_flip() + ylab("Type 1 Error") + xlab("Sufficiency Threshold") +
  scale_y_continuous(limits=c(0, 0.5), breaks=c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)) +
  scale_shape_manual(values=c(1,2,15,16,17))+
  scale_color_manual(values=c('#999999','#E69F00', '#56B4E9', '#D9717D','#BECA55')) +
  geom_hline(yintercept=c(0.05, 0.1),color="gray",alpha=0.7) + 
  theme_minimal() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        legend.position = "right", legend.justification = "center",
        legend.title = element_text(size=5),
        legend.text = element_text(size=5),
        legend.key.size = unit(0.25, "cm"),
        legend.key.width = unit(0.5,"cm"),
        axis.text=element_text(size=4),
        axis.title=element_text(size=5),
        strip.text.x = element_text(size = 4),
        strip.text.y = element_text(size = 4)) +
  facet_wrap(~SampleSize,ncol=1)

Type1err_graph

#combine graphs 
Combined_graph <- ggarrange(Type1err_graph, Power_graph, labels = c("A", "B"), ncol = 2, nrow = 1,common.legend=T)
Combined_graph
ggsave("SampleSize_Type1ErrPower_Prior.png")


