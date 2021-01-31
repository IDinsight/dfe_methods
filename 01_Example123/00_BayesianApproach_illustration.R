
###########################
### 00_BayesianApproach_illustration
###########################

# Date: 10 DEC 2020 
# Author: Torben Fischer 
# Contact: torben.fischer@idinsight.org
# Edits: 

# Description: 
#   This routine illustrate Bayesian Priors and Posteriors

setwd("C:\Users\tofis\Dropbox (IDinsight)\IE Design for DFEs\Code\03_Output")
library(tidyr)


quantiles <- seq(0, 1, length.out = 100)  # Range of probability / success parameter
set.seed(22)
# for each success probability calculate likeliness
likel <- data.frame(quantiles,dbinom(x=3, prob = quantiles, size = 10) ,"Likelihood")
  colnames(likel) <- c("Range","Density","Distribution")
# uninformative
set.seed(22)
prior_uninform <- data.frame(quantiles,dnorm(x=quantiles,mean=0.5,sd=1.5),"Prior" )
colnames(prior_uninform) <- c("Range","Density","Distribution")
  
post_unin_ns = likel[2]*prior_uninform[2]  
post_unin_s = post_unin_ns/sum(post_unin_ns)
post_uninform_noscale=data.frame(quantiles,post_unin_ns,"Posterior (noscale)")
  colnames(post_uninform_noscale) <- c("Range","Density","Distribution")
post_uninform=data.frame(quantiles,post_unin_s,"Posterior")
  colnames(post_uninform) <- c("Range","Density","Distribution")  

sum(post_uninform_noscale[2])
sum(post_uninform[2]) 
  
# informative
set.seed(22)
prior_inform <- data.frame(quantiles,dnorm(x=quantiles,mean=0.5,sd=0.1)/10,"Informative Prior")
  colnames(prior_inform) <- c("Range","Density","Distribution")  

post_in_ns = likel[2]*prior_inform[2]  
post_in_s = post_in_ns/sum(post_in_ns)  
post_inform_noscale=data.frame(quantiles,post_in_ns,"Posterior (noscale)")
  colnames(post_inform_noscale) <- c("Range","Density","Distribution")
post_inform=data.frame(quantiles,post_in_s,"Posterior")  
  colnames(post_inform) <- c("Range","Density","Distribution")

sum(post_inform_noscale[2])
sum(post_inform[2])


  
# plot data   
df_inform <- rbind(likel,prior_inform,post_inform_noscale,post_inform)
df_uninform <- rbind(likel,prior_uninform,post_uninform_noscale,post_uninform)

# informative plot
summary(df_inform)
inform_gr <- ggplot(data=df_inform, aes(x=Range,y=Density,group=factor(Distribution),color=Distribution )) + 
  geom_line() + xlab("P(Head)") +
  scale_y_continuous(limits=c(0, 0.5), breaks=seq(0,0.5,0.05)) +
  theme_bw()

inform_gr
ggsave(file.path(path(getwd(),"BayesInformative.png",ext="")), inform_gr, height=7, width=7)

# uninformative plot
summary(df_uninform)
uninform_gr <- ggplot(data=df_uninform, aes(x=Range,y=Density,group=factor(Distribution),color=Distribution )) + 
  geom_line() + xlab("P(Head)") +
  scale_y_continuous(limits=c(0, .5), breaks=seq(0,.5,0.05)) +
  theme_bw()

uninform_gr
ggsave(file.path(path(getwd(),"BayesUninformative.png",ext="")), uninform_gr, height=7, width=7)

