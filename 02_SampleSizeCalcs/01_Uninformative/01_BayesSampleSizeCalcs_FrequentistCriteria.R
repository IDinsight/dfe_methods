


#####################################
### Bayesian Sample Size Calculations - Using Frequentist Operating Characteristics 
#####################################


# Date: 18 DEC 2019 
# Author: Torben Fischer 
# Contact: torben.fischer@idinsight.org
# Edits: 28 JAN 2020 - added several priors, added double-barreled t-test 

# Description: 
#   This routine calculates the required sample sizes for a Bayesian approach using frequentist operating criteria
#   in a comaring to a benchmark effect size scenario.
#   The frequentist criteria used to assess the quality of the study are: 
#   1) (Freq) Type 1 error: proportion of simulation studies that suggest to take the "wrong" decision when assuming 
#      an effect size of 0. 
#   2) (Freq) Statistical Power: proportion of simulaton studies that suggest to take the "right" decision when assuming 
#      the smallest effect size for which that decision is correct. 

# Content: 
#   0. Set-up
#   1. Bayesian Sample Size Calculations 
#       i) Specifying design & analysis priors, parameter assumptions 
#      ii) Defining a decision rule 
#     iii) Data Generating process 
#      iv) Frequentist Type 1 Error Rate 
#       v) Frequentist Statistical Power  

##############################################################


#   0. Set-up
##############################################################

rm(list=ls()) # clear current environment

# Working directory
setwd("C:/Users/tofis/Dropbox (IDinsight)/IE Design for DFEs/Code") # working directory, specify your path here. 
library(fs) # use of path()
outpath = path(getwd(),"03_Output","01_SampleSizeCalcs",ext="") # define output path folder,make sure sub-folder exists
output_name = "BayesSampleSizeCalcs_freqCritOpsstat_N5000only.RData"
# Elaborate folder structure if special 

# Load packages 

library(rstan)
# Stan Options
rstan_options(auto_write = TRUE)             # avoid recompilation of models
options(mc.cores = parallel::detectCores()-1)  # parallelize across all CPUs
Sys.setenv(LOCAL_CPPFLAGS = '-march=native') # improve execution time

library(rstudioapi) # to fit rstan

library(future.apply) # parellize simulations

#install.packages("tidyverse")    
library("tidyverse")    ## load the core tidyverse packages, incl. dplyr
library(dplyr)
library(magrittr) 
library(tidyr)
library(wesanderson)
library("ggpubr") ## combine graphs


#   1. Bayesian Sample Size Calculations 
##############################################################


#       i) Specifying priors & assumptions
##########################################

# Analysis Priors: 
# b ~ N(.1,.1^2)
# a ~ N(-1.7,2^2)

# Assumptions:

# i) assume fixed variance sigma for the likelihood:
sigma = 1 
#   --> the standard deviation of height-for-age is generally around 1 for most populations 
#   --> thus this is a fairly benign assumption.

# ii) b = 0 for Type 1 Error Calculation 
b_t1e = 0

# iii) b = cash benchmark + eps for Power calculation 
bnmk = 0.2   # benchmark effect size to compare to
EPS = 0.05    # deviation from benchmark, tolerance for Type 1 error calc  

b_power = bnmk + EPS 

# iv) Simulation sample sizes
RS = 1000     # Type 1 error calcl No. simulation samples to draw for a given sample size
REPS = 1000    # No. of sims for Statistical Power Calcs, set to 1000 

#      ii) Defining a decision rule 
###################################

# Rule: Support/fund/scale program if P(ATE>bnmk)>cert_thresh
# In words: if more than "cert_thresh*100%" of posterior probability mass is on values larger than "bnmk", then support. 
# Other more sophisticated decision rules exist, see text. 

cert_thresh = c(0.5,0.6,0.7,0.8) # vector of certainty thresholds, probability to be passed to take decision 

#     iii) Data Generating process 
###################################

fake_data = function(input){  
  SEED = input[1]
  a = input[2]
  b = input[3]
  SIGMA = input[4]  
  n = input[5]
  
  # n - sample size
  # a - constant (drawn from prior) 
  # b - ATE (drawn from prior or fixed)
  set.seed(SEED)
  eps = rnorm(n,mean=0,sd=SIGMA)
  set.seed(SEED)
  T = rbinom(n,1,0.5) # random draw of treatment status (binary)
  y = a + b*T + eps   # outcome data, assume simple linear and separable 
  fak_dat = data.frame(y,T) 
  return(fak_dat)
}  


#      iv) Type 1 Error Rate 
###################################

# set-up data simulation process
#s = c(50,100,200,300,500,1000,2000,4000,5000) # array of sample sizes to consider
s = c(5000) # array of sample sizes to consider
seed = 20200128       # initial seed to draw from DGP 

seed.setter = function(startseed,reps){ 
  # Purpose: Generates seeds for different draws based on initial seed
  # startseed = initial seed
  # reps = # of different seeds / samples required
  seeds = seed + c((0:(reps-1)))
  return(seeds)
}

# Generate draws from the (design) prior distribution for parameters other than b 
#     (b is fixed, note: sigma also fixed, see step i)
set.seed(seed)
a_priordraws = rnorm(length(s)*RS,mean=-1.7,sd=2)
length(a_priordraws)

# Generate seeds to draw from DGP 
SEEDS = as.list(seed.setter(seed,RS*length(s)))
length(SEEDS)

# Combine moving input pieces
INPUTS = cbind( unlist(SEEDS),a_priordraws,rep(b_t1e,length(s)*RS),rep(sigma,length(s)*RS),sort(rep(s,RS)) ) 
  # INPUTS must have same format as input function

INPUTS = lapply(seq_len(nrow(INPUTS)), function(i) INPUTS[i,]) # turn into list argument 

# Draw from DGP 
dat = list()  # data for all sample sizes 
dat = lapply(FUN=fake_data,X=INPUTS)

do.call(rbind,lapply(FUN=dim,X=dat)) # check dimensions

#dat1 = lapply(FUN=fake_data,X=SEEDS,a=-1.4,b=bnmk,tol=eps,n=s[1])  # This is the data for one sample size 

# Specify Stan model (as separate file)
benchmarking_uninformative.mod <- stan_model(file = path(getwd(),"02_sampleSizeCalcs","01_Uninformative","01_benchmarking_uninformative.stan"))   ## C++ Compilation


################# Test: Estimate Bayesian Model for a given data set, get decision for that data set 
#####################################################################################################

# Specify Data (this will be done through parallized apply function later)
stan.dat = list(y=dat[[1]][,1],
                treat=dat[[1]][,2],
                N=length(dat[[1]][,1])
)  
# Note this is a list of lists (varies over sample size and # draws from DGP)

# Model Fitting 
lm.est <- sampling(benchmarking_uninformative.mod,    # compiled model
                   data = stan.dat,                   # data input
                   algorithm = "NUTS",                # algorithm
                   control = list(                    # control arguments
                     adapt_delta = .85
                   ),
                   save_warmup = FALSE,               # discard warmup sims
                   sample_file = NULL,                # no sample file --> usually recommended to use
                   diagnostic_file = NULL,            # no diagnostic file --> usually recommended to use 
                   pars = c("beta", "alpha"),         # select parameters
                   iter = 2000L,                      # iter per chain
                   warmup = 1000L,                    # warmup period
                   thin = 2L,                         # thinning factor
                   chains = 4L,                       # num. chains
                   cores = 4L,                        # num. cores
                   seed = seed)                       # seed

# str(lm.est) # model summary
# Model descision 

print(lm.est) # summarize model 
extracted_values<-rstan::extract(lm.est, permuted = TRUE) # extract the posterior draws from the fit object

# estimate the probability that beta is greater than bnmk
prob_better <- sum(extracted_values$beta >bnmk)/length(extracted_values$beta) # bnmk is 'external'
decision = prob_better >= cert_thresh # cert_thresh is 'external'

# Frequentist decision rule 
# Test 1 
test1 = t.test(x=dat[[1]][dat[[1]]$T==1,1],y=dat[[1]][dat[[1]]$T==0,1],mu=0,alternative="greater",conf.level = 0.95,data=dat[[1]])
freq_t1_res = (test1$p.value < (1 - attr(test1$conf.int,"conf.level")))
# Test 2 
test2 = t.test(x=dat[[1]][dat[[1]]$T==1,1],y=dat[[1]][dat[[1]]$T==0,1],mu=0.2,alternative="greater",conf.level = 0.8,data=dat[[1]])
freq_t2_res = (test2$p.value < (1 - attr(test2$conf.int,"conf.level")))

######################### END - Test for 1 dataset ######################

#### Investigate posterior estimate for beta
x= extracted_values$beta
h<-hist(x, breaks=10, col="gray", xlab="Beta",
        main="Posterior Density against Normal Curve") ### N=50
## add a normal: 
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=2) 
d = density(x) # kernel density
d$Y = d$y*diff(h$mids[1:2]*length(x))
lines(d$x,d$Y,col="green",lwd=2)




# Function that takes a dataset, pre-defined Stan model, decision rule 
#       --> output: decision for a given 

decision_func = function(stanMod,bnmk,cert_thresh,fitData){
  # Description: 
  # Given a model, data set and decision rule, this function fits the model and determines whether a decision 
  # should be taken or not. 
  # Inputs: 
  # stanMod: Rstan model/file to be fitted 
  # fitData: data to be used for estimation, same structure as specified in stanMod 
    # here: fitData should be a data frame with 2 columns (y,T)
  # bnmk: benchmark effect size, specified above
  # cert_thresh: certainty threshold to take decision, specified above. 
  # Output: 
  # prob_better: probability that intervention is better than decision. 
  # decision: binary indicator for whether benchmark comparison is successful
  
  stan.dat = list(y=fitData[,1],
                  treat=fitData[,2],
                  N=length(fitData[,1])
  )  
  # Note this is a list of lists (varies over sample size and # draws from DGP)
  
  # Fit Model 
  lm.est <- sampling(stanMod,                           # compiled model
                     data = stan.dat,                   # data input
                     algorithm = "NUTS",                # algorithm
                     control = list(                    # control arguments
                       adapt_delta = .85
                     ),
                     save_warmup = FALSE,               # discard warmup sims
                     sample_file = NULL,                # no sample file --> usually recommended to use
                     diagnostic_file = NULL,            # no diagnostic file --> usually recommended to use 
                     pars = c("beta", "alpha"),         # select parameters
                     iter = 2000L,                      # iter per chain
                     warmup = 1000L,                    # warmup period
                     thin = 2L,                         # thinning factor
                     chains = 4L,                       # num. chains
                     cores = 4L,                        # num. cores
                     seed = seed)                       # seed
  
  # Draw from posterior 
  extracted_values<-rstan::extract(lm.est, permuted = TRUE) # extract the posterior draws from the fit object
  
  # estimate the probability that ATE > bnmk as % of draws from posterior
  prob_better <- sum(extracted_values$beta >bnmk)/length(extracted_values$beta) # bnmk is 'external'
  #decision = prob_better >= cert_thresh # cert_thresh is 'external
  #Note: can be inferred post-analysis based on prob-better
  
  # Frequentist decision rule 
  # Test 1 
  test1 = t.test(x=fitData[fitData$T==1,1],y=fitData[fitData$T==0,1],mu=0,alternative="greater",conf.level = 0.95,data=fitData)
  freq_t1_res = (test1$p.value < (1 - attr(test1$conf.int,"conf.level")))
  # Test 2 
  test2 = t.test(x=fitData[fitData$T==1,1],y=fitData[fitData$T==0,1],mu=0.2,alternative="greater",conf.level = 0.8,data=fitData)
  freq_t2_res = (test2$p.value < (1 - attr(test2$conf.int,"conf.level")))
  
  # Output
  out = cbind(prob_better,freq_t1_res,freq_t2_res,stan.dat[[3]])
  colnames(out)=c("prob_betbnmk","Reject Test1","Reject Test2","samp_size")
  return(out)
}


#decision_func(stanMod=benchmarking_uninformative.mod,bnmk=bnmk,cert_thresh=cert_thresh,fitData=dat[[1]][1][[1]])

# Run the model for all samples of a given sample size
start_time_np <- Sys.time()

res=lapply(X=dat,FUN=decision_func,stanMod=benchmarking_uninformative.mod,bnmk=bnmk,cert_thresh=cert_thresh) 
decisions=do.call(rbind,res)

end_time_np <- Sys.time()
end_time_np - start_time_np # 6 sample sizes with 10 samples each take ~5minutes. 

dec = as.data.frame(decisions)
dec = cbind(dec,dec[,2]*dec[,3])
names(dec)=c("prob_betbnmk","Freq_test1", "Freq_test2","samp_size","Freq_dec")

####################################
# iv) Frequentist Type 1 Error Rate
####################################

cert_thresh

type1err = dec %>% 
  group_by(samp_size) %>%
  summarize(type1err_50 = mean(prob_betbnmk>cert_thresh[1], na.rm=TRUE),
            type1err_60 = mean(prob_betbnmk>cert_thresh[2], na.rm=TRUE),
            type1err_70 = mean(prob_betbnmk>cert_thresh[3], na.rm=TRUE),
            type1err_80 = mean(prob_betbnmk>cert_thresh[4], na.rm=TRUE),
            freq_barrel = mean(Freq_dec,na.rm=TRUE),
            freq_test1 = mean(Freq_test1,na.rm=TRUE),
            freq_test2 = mean(Freq_test2,na.rm=TRUE)) 
  

# Select sample sizes of interest (to test for power)
selected_sams = type1err[type1err[,2]<=0.05,1]
Sams = dim(selected_sams)[1]


## Graphical Illustration
names(type1err)
graph_type1err <- type1err[,1:5] # keep Bayesian cols 
colnames(graph_type1err)=c("Sample Size","50","60","70","80")
graph_type1err_long <- pivot_longer(graph_type1err,cols=2:5,names_to="Threshold",values_to= "Type 1 Error")
graph_type1err_long = as.data.frame(graph_type1err_long)
rep2 = as.numeric(graph_type1err_long[,2])
graph_type1err_long=as.data.frame(cbind(graph_type1err_long[,1],rep2,graph_type1err_long[,3]))
colnames(graph_type1err_long)=c("SampSize","Threshold","Type1Error")

# Mannually add estimates for SampeSize 5000
# # A tibble: 1 x 8
# samp_size type1err_50 type1err_60 type1err_70 type1err_80 freq_barrel freq_test1 freq_test2
# <dbl>       <dbl>       <dbl>       <dbl>       <dbl>       <dbl>      <dbl>      <dbl>
#   1      5000           0           0           0           0           0       0.05          0

n5000=as.data.frame(matrix(c(rep(5000,4),50,60,70,80,0,0,0,0),nrow=4,ncol=3,byrow=F))
colnames(n5000)=c("SampSize","Threshold","Type1Error")

graph_type1err_long=rbind(graph_type1err_long,n5000)


# Set axis limits c(min, max)
min <- 0
max <- 0.28
ncols=length(unique(graph_type1err_long$SampSize))

Type1Err_graph <- graph_type1err_long %>%
  ggplot( aes(x=Threshold, y=Type1Error, group=factor(SampSize), color=factor(SampSize)))  +
  geom_line(size=1) + geom_point(size=2) +  
  ylab("Type 1 Error") + xlab("Sufficiency Threshold (%)") +
  scale_y_continuous(limits=c(min, max), breaks=c(0, 0.05, 0.1, 0.15, 0.2, 0.25)) +
  geom_hline(yintercept = 0.05) +   scale_colour_discrete("Sample Size") +
  geom_hline(yintercept = seq(min,max,0.025),linetype=3, size=0.5) +
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), 
          legend.position = "right", legend.justification = "center",
          legend.title = element_text(size=8),
          legend.text = element_text(size=6),
          legend.key.size = unit(0.3, "cm"),
          legend.key.width = unit(0.5,"cm"),
          axis.text=element_text(size=6),
          axis.title=element_text(size=8))

Type1Err_graph

###################################
# v) Frequentist Statistical Power 
###################################

if(Sams==0){
  error1 = "None of the sample sizes specified has favourable type 1 error rates. Select other."
  print(error1)
}

Sams=8
#selected_sams = c(50,100,200,300,500,1000,2000,4000,5000) # For illustrative purposes, can be integrated with sample size from T1E calculation
selected_sams = c(5000) 

# step 1 - set-up
seed = 19901224       # initial seed to draw from DGP 

# Generate draws from the (design) prior distribution for parameters other than b (as fixed, note: sigma also fixed, see step i)
set.seed(seed)
a_priordraws = rnorm(length(selected_sams)*REPS,mean=-1.7,sd=2)
length(a_priordraws)

# Generate seeds to draw from DGP 
SEEDS = as.list(seed.setter(seed,REPS*length(selected_sams)))
length(SEEDS)

# Combine moving input pieces
powerINPUTS = cbind( unlist(SEEDS),a_priordraws,rep(b_power,length(selected_sams)*REPS),
                rep(sigma,length(selected_sams)*REPS),sort(rep(selected_sams,REPS)) ) 
# INPUTS must have same format as input function

powerINPUTS = lapply(seq_len(nrow(powerINPUTS)), function(i) powerINPUTS[i,]) # turn into list argument 

# Draw from DGP 
powerdat = list()  # data for all sample sizes 
powerdat = lapply(FUN=fake_data,X=powerINPUTS)

do.call(rbind,lapply(FUN=dim,X=powerdat)) # check dimensions


# Simulate studies
start_time_np <- Sys.time()

power_res=lapply(X=powerdat,FUN=decision_func,stanMod=benchmarking_uninformative.mod,bnmk=bnmk,cert_thresh=cert_thresh) 
power_decisions=do.call(rbind,power_res)

end_time_np <- Sys.time()
end_time_np - start_time_np # 6 sample sizes with 10 samples each take ~5minutes. 

power_dec = as.data.frame(power_decisions)
power_dec = cbind(power_dec,power_dec[,2]*power_dec[,3])
names(power_dec)=c("prob_betbnmk","Freq_test1", "Freq_test2","samp_size","Freq_dec")

##########################################################
# calculate frequentist statistical power (by sample size)
##########################################################

cert_thresh

freq_power = power_dec %>% 
  group_by(samp_size) %>%
  summarize(power_50 = mean(prob_betbnmk>cert_thresh[1], na.rm=TRUE),
            power_60 = mean(prob_betbnmk>cert_thresh[2], na.rm=TRUE),
            power_70 = mean(prob_betbnmk>cert_thresh[3], na.rm=TRUE),
            power_80 = mean(prob_betbnmk>cert_thresh[4], na.rm=TRUE),
            freq_power_barrel = mean(Freq_dec,na.rm=TRUE),
            freq_power_test1 = mean(Freq_test1,na.rm=TRUE),
            freq_power_test2 = mean(Freq_test2,na.rm=TRUE)) 

## Graphical Illustration
names(freq_power)
graph_freq_power <- freq_power[,1:5] # keep Bayesian cols 
colnames(graph_freq_power)=c("Sample Size","50","60","70","80")
graph_freq_power_long <- pivot_longer(graph_freq_power,cols=2:5,names_to="Threshold",values_to= "Power")
graph_freq_power_long = as.data.frame(graph_freq_power_long)
rep2 = as.numeric(graph_freq_power_long[,2])
graph_freq_power_long=as.data.frame(cbind(graph_freq_power_long[,1],rep2,graph_freq_power_long[,3]))
colnames(graph_freq_power_long)=c("SampSize","Threshold","Power")

# manual add the results for n=5000
# A tibble: 1 x 8
# samp_size power_50 power_60 power_70 power_80 freq_power_barrel freq_power_test1 freq_power_test2
# <dbl>    <dbl>    <dbl>    <dbl>    <dbl>             <dbl>            <dbl>            <dbl>
#  1      5000    0.962    0.934    0.895    0.819             0.823                1            0.823

n5000=as.data.frame(matrix(c(rep(5000,4),50,60,70,80,0.962,0.934,0.895,0.819),nrow=4,ncol=3,byrow=F))
colnames(n5000)=c("SampSize","Threshold","Power")

graph_freq_power_long=rbind(graph_freq_power_long,n5000)


# Set axis limits c(min, max)
min <- 0.24
max <- 1
ncols=length(unique(graph_freq_power_long$SampSize))

Power_graph <- graph_freq_power_long %>%
  ggplot( aes(x=Threshold, y=Power, group=factor(SampSize), color=factor(SampSize)))  +
  geom_line(size=1) + geom_point(size=2) +  
  ylab("Statistical Power") + xlab("Sufficiency Threshold (%)") +
  scale_y_continuous(limits=c(min, max), breaks=c(0.3, 0.4, 0.5, 0.6, 0.7,0.8,0.9,1)) +
  geom_hline(yintercept = 0.8) +   scale_colour_discrete("Sample Size") +
  geom_hline(yintercept = seq(0.3,max,0.1),linetype=3, size=0.5) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        legend.position = "right", legend.justification = "center",
        legend.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.3, "cm"),
        legend.key.width = unit(0.5,"cm"),
        axis.text=element_text(size=6),
        axis.title=element_text(size=8))

Power_graph

# combine graphs
Combined_graph <- ggarrange(Type1Err_graph, Power_graph, labels = c("A", "B"), ncol = 2, nrow = 1, align="hv",common.legend=T)
Combined_graph
ggsave(file.path(outpath,"SampleSize_Type1ErrPower.png"), Combined_graph)

# save and exit 
save.image(file=path(outpath,"00_Rdata",output_name))

write.csv(type1err,file=path(outpath,"02_Type1Err","type1err.csv",ext=""))
write.csv(freq_power,file=path(outpath,"01_Power","power.csv",ext=""))


##### Drawing normal posterior and test distribution 
####################################################

n=10000
qnorm(.2,mean=0,sd=1)

pnorm(-0.8416212,mean=0,sd=1)








