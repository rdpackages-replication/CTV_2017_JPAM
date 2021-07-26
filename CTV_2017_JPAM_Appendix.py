################################################################################
# Comparing Inference Approaches for RD Designs:
# A Reexamination of the Effect of Head Start on Child Mortality
# Supplemental Appendix
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
# Python code created by Ricardo Masini
# Last update: 26-JUL-2021
################################################################################
## SOFTWARE WEBSITE: https://rdpackages.github.io/
################################################################################
## TO INSTALL/DOWNLOAD PYTHON PACKAGE:
## RDROBUST: pip install rdrobust
################################################################################
## NOTE: if you are using RDROBUST version 2020 or newer, the option 
## masspoints="off" and stdvars=True may be needed to replicate the results
## in the paper.
################################################################################
## NOTE: if you are using RDDENSITY version 2020 or newer, the option 
## masspoints=False may be needed to replicate the results in the paper.
################################################################################

from rdrobust import rdrobust, rdplot
import numpy as np
import pandas as pd

#############################################################
## Load and setup data
#############################################################

data = pd.read_csv('headstart.csv')

# 1960 Census covariates

cutoff = 59.1984

X60 = data[['census1960_pop',
            'census1960_pctsch1417',
            'census1960_pctsch534', 
            'census1960_pctsch25plus',
            'census1960_pop1417',
            'census1960_pop534',
            'census1960_pop25plus',
            'census1960_pcturban',
            'census1960_pctblack']]

# 1990 Census covariates

X90 = data[['census1990_pop',
            'census1990_pop1824',
            'census1990_pop2534',
            'census1990_pop3554',
            'census1990_pop55plus',
            'census1990_pcturban',
            'census1990_pctblack',
            'census1990_percapinc']]

# Outcome, normalized running variable and treatment

Y = data['mort_age59_related_postHS']
Rraw = data['povrate60']
R = data['povrate60'] - cutoff

D = ((R>=0).values)*1

# Placebo outcomes

Plac = data[['mort_age59_injury_postHS',
             'mort_age59_all_postHS',
             'mort_age25plus_related_postHS',
             'mort_age25plus_injuries_postHS', 
             'mort_age59_related_preHS',
             'mort_wh_age59_related_postHS',
             'mort_bl_age59_related_postHS']]

###################################################################
## Table SA-1: Outcome and Score Descriptive Stats
###################################################################

# TableSA1 <- array(NA,dim=c(10,2))

# RY <- cbind(Rraw,Y)

# TableSA1[1,] <- colMeans(RY,na.rm=T)
# TableSA1[2,] <- apply(RY,2,sd,na.rm=T)
# TableSA1[3,] <- apply(RY,2,min,na.rm=T)
# TableSA1[4:8,] <- apply(RY,2,function(x) quantile(x,probs=c(0.1,.25,0.5,.75,0.9),na.rm=T))
# TableSA1[9,] <- apply(RY,2,max,na.rm=T)
# TableSA1[10,] <- apply(RY,2,function(x) sum(!is.na(x)))

# round(TableSA1,3)

###################################################################
## Figure SA-1 and SA-2: Histograms
###################################################################

# hist(Rraw,breaks=60)
# hist(Y,breaks=60)
# hist(Y[Y>0],breaks=60)

# w <- 1
# ww <- which(abs(R)<=w)
# hist(Rraw[ww],breaks=20)
# hist(Y[ww],breaks=20)
# hist(Y[ww][Y[ww]>0],breaks=12)

# w <- 3
# ww <- which(abs(R)<=w)
# hist(Rraw[ww],breaks=50)
# hist(Y[ww],breaks=60)
# hist(Y[ww][Y[ww]>0],breaks=60)

# w <- 9
# ww <- which(abs(R)<=w)
# hist(Rraw[ww],breaks=50)
# hist(Y[ww],breaks=60)
# hist(Y[ww][Y[ww]>0],breaks=60)


###################################################################
## Figure SA-3: Rdplots
###################################################################

rdplot(Y,Rraw,c=cutoff,y_lim=(0,9))

for k in Plac.columns: 
  rdplot(Plac[k],Rraw,c=cutoff,y_lim=(np.min(Y),np.max(Y)))


###################################################################
## Tables SA-2 to SA-5: Summary statistics and Difference in means
###################################################################

# TableSA2 <- array(NA,dim=c(18,9))
# TableSA3 <- array(NA,dim=c(18,9))
# TableSA4 <- array(NA,dim=c(18,9))
# TableSA5 <- array(NA,dim=c(18,9))

# alldata <- cbind(Y,Rraw,Plac,X60)

# # SA-2

# w <- 100

# for(col in 1:ncol(alldata)){
  
#   var <- alldata[,col]
#   w.t <- which(abs(R)<=w & complete.cases(var) & D==1)
#   w.c <- which(abs(R)<=w & complete.cases(var) & D==0)
#   w.a <- which(abs(R)<=w & complete.cases(var))
  
#   n.t <- length(var[w.t])
#   m.t <- mean(var[w.t])
#   se.t <- sd(var[w.t])/sqrt(length(var[w.t]))
  
#   n.c <- length(var[w.c])
#   m.c <- mean(var[w.c])
#   se.c <- sd(var[w.c])/sqrt(length(var[w.c]))
  
#   t.stat <- abs(m.t-m.c)/sqrt(se.t^2+se.c^2)
  
#   p.val <- 2*(1-pnorm(t.stat))
  
#   TableSA2[col,1] <- n.c
#   TableSA2[col,2] <- m.c
#   TableSA2[col,3] <- se.c
#   TableSA2[col,4] <- n.t
#   TableSA2[col,5] <- m.t
#   TableSA2[col,6] <- se.t
#   TableSA2[col,7] <- m.t-m.c
#   TableSA2[col,8] <- sqrt(se.t^2+se.c^2)
#   TableSA2[col,9] <- p.val
  
# }
# round(TableSA2,3)

# # SA-3

# w <- 9

# for(col in 1:ncol(alldata)){
  
#   var <- alldata[,col]
#   w.t <- which(abs(R)<=w & complete.cases(var) & D==1)
#   w.c <- which(abs(R)<=w & complete.cases(var) & D==0)
#   w.a <- which(abs(R)<=w & complete.cases(var))
  
#   n.t <- length(var[w.t])
#   m.t <- mean(var[w.t])
#   se.t <- sd(var[w.t])/sqrt(length(var[w.t]))
  
#   n.c <- length(var[w.c])
#   m.c <- mean(var[w.c])
#   se.c <- sd(var[w.c])/sqrt(length(var[w.c]))
  
#   t.stat <- abs(m.t-m.c)/sqrt(se.t^2+se.c^2)
  
#   p.val <- 2*(1-pnorm(t.stat))
  
#   TableSA3[col,1] <- n.c
#   TableSA3[col,2] <- m.c
#   TableSA3[col,3] <- se.c
#   TableSA3[col,4] <- n.t
#   TableSA3[col,5] <- m.t
#   TableSA3[col,6] <- se.t
#   TableSA3[col,7] <- m.t-m.c
#   TableSA3[col,8] <- sqrt(se.t^2+se.c^2)
#   TableSA3[col,9] <- p.val
  
# }
# round(TableSA3,3)

# SA-4

# w <- 3

# for(col in 1:ncol(alldata)){
  
#   var <- alldata[,col]
#   w.t <- which(abs(R)<=w & complete.cases(var) & D==1)
#   w.c <- which(abs(R)<=w & complete.cases(var) & D==0)
#   w.a <- which(abs(R)<=w & complete.cases(var))
  
#   n.t <- length(var[w.t])
#   m.t <- mean(var[w.t])
#   se.t <- sd(var[w.t])/sqrt(length(var[w.t]))
  
#   n.c <- length(var[w.c])
#   m.c <- mean(var[w.c])
#   se.c <- sd(var[w.c])/sqrt(length(var[w.c]))
  
#   t.stat <- abs(m.t-m.c)/sqrt(se.t^2+se.c^2)
  
#   p.val <- 2*(1-pnorm(t.stat))
  
#   TableSA4[col,1] <- n.c
#   TableSA4[col,2] <- m.c
#   TableSA4[col,3] <- se.c
#   TableSA4[col,4] <- n.t
#   TableSA4[col,5] <- m.t
#   TableSA4[col,6] <- se.t
#   TableSA4[col,7] <- m.t-m.c
#   TableSA4[col,8] <- sqrt(se.t^2+se.c^2)
#   TableSA4[col,9] <- p.val
  
# }
# round(TableSA4,3)

# SA-5

# w <- 1

# for(col in 1:ncol(alldata)){
  
#   var <- alldata[,col]
#   w.t <- which(abs(R)<=w & complete.cases(var) & D==1)
#   w.c <- which(abs(R)<=w & complete.cases(var) & D==0)
#   w.a <- which(abs(R)<=w & complete.cases(var))
  
#   n.t <- length(var[w.t])
#   m.t <- mean(var[w.t])
#   se.t <- sd(var[w.t])/sqrt(length(var[w.t]))
  
#   n.c <- length(var[w.c])
#   m.c <- mean(var[w.c])
#   se.c <- sd(var[w.c])/sqrt(length(var[w.c]))
  
#   t.stat <- abs(m.t-m.c)/sqrt(se.t^2+se.c^2)
  
#   p.val <- 2*(1-pnorm(t.stat))
  
#   TableSA5[col,1] <- n.c
#   TableSA5[col,2] <- m.c
#   TableSA5[col,3] <- se.c
#   TableSA5[col,4] <- n.t
#   TableSA5[col,5] <- m.t
#   TableSA5[col,6] <- se.t
#   TableSA5[col,7] <- m.t-m.c
#   TableSA5[col,8] <- sqrt(se.t^2+se.c^2)
#   TableSA5[col,9] <- p.val
  
# }
# round(TableSA5,3)

###################################################################
## Figure SA-4: Falsification Test - Rdplots
###################################################################

rdplot(Y,Rraw,c=cutoff,binselect='es',y_lim=(0,8))
rdplot(Y,Rraw,c=cutoff,binselect='qs',y_lim=(0,8))
rdplot(Y,Rraw,c=cutoff,binselect='qsmv',y_lim=(0,8))

###################################################################
## Tables SA-6 to SA-8: Local Polynomial, Main and Placebo Outcomes
###################################################################

# Table SA-6
mataux = pd.concat([Y, Plac], axis=1)

# Table SA-6.0:
for col in mataux.columns:# Table SA-6.0:
    rdrobust(mataux[col],R,p=0,q=1,bwselect='cerrd')
    rdrobust(mataux[col],R,p=0,q=1,bwselect='mserd')
    rdrobust(mataux[col],R,p=0,q=1,h=9)
    rdrobust(mataux[col],R,p=0,q=1,h=18)
    
# Table SA-6.1:
for col in mataux.columns:
    rdrobust(mataux[col],R,p=1,bwselect='cerrd')
    rdrobust(mataux[col],R,p=1,bwselect='mserd')
    rdrobust(mataux[col],R,p=1,h=9)
    rdrobust(mataux[col],R,p=1,h=18)


# Table SA-6.2:
rdrobust(Y,R,p=0,q=1,bwselect='cerrd')
rdrobust(Y,R,p=0,q=1,bwselect='mserd')
rdrobust(Y,R,p=0,q=1,h=9)
rdrobust(Y,R,p=0,q=1,h=18)
for col in Plac.columns:
    rdrobust(Plac[col],R,p=0,q=1,bwselect='cerrd')
    rdrobust(Plac[col],R,p=0,q=1,bwselect='mserd')
    rdrobust(Plac[col],R,p=0,q=1,h=9)
    rdrobust(Plac[col],R,p=0,q=1,h=18)

# Table SA-6.3:
rdrobust(Y,R,p=1,bwselect='cerrd')
rdrobust(Y,R,p=1,bwselect='mserd')
rdrobust(Y,R,p=1,h=9)
rdrobust(Y,R,p=1,h=18)
for col in Plac.columns:
    rdrobust(Plac[col],R,p=1,bwselect='cerrd')
    rdrobust(Plac[col],R,p=1,bwselect='mserd')
    rdrobust(Plac[col],R,p=1,h=9)
    rdrobust(Plac[col],R,p=1,h=18)

# Table SA-6.4:
rdrobust(Y,R,p=4,q=5,bwselect='cerrd')
rdrobust(Y,R,p=4,q=5,bwselect='mserd')
rdrobust(Y,R,p=4,q=5,h=9)
rdrobust(Y,R,p=4,q=5,h=18)

for col in Plac.columns:
    rdrobust(Plac[col],R,p=4,q=5,bwselect='cerrd')
    rdrobust(Plac[col],R,p=4,q=5,bwselect='mserd')
    rdrobust(Plac[col],R,p=4,q=5,h=9)
    rdrobust(Plac[col],R,p=4,q=5,h=18)

###################################################################
## Figure SA-5: Local Randomization Methods - Window Selection
###################################################################

# wreps <- 1000

# ## NOTE: the plots are drawn using the asymptotic p-value to speed up the process.
# ## Remove the "approx" option to use randinf and replicate the results in the paper.

# tmp <- rdwinselect(R,cbind(Plac[,1],X60),reps=wreps,stat='ksmirnov',wmin=.3,wstep=.2,nwin=40,level=.2,quietly=TRUE,approx=TRUE)
# plot(tmp$results[,7],tmp$results[,1],ylab='p-values',xlab='bandwidth')
# abline(v=1.1,lty='dashed')

# tmp <- rdwinselect(R,cbind(Plac[,1],X60),reps=wreps,stat='ttest',wmin=.3,wstep=.2,nwin=40,level=.2,quietly=TRUE,approx=TRUE)
# plot(tmp$results[,7],tmp$results[,1],ylab='p-values',xlab='bandwidth')
# abline(v=1.5,lty='dashed')

# tmp <- rdwinselect(R,cbind(Plac[,1],X60),reps=wreps,stat='ranksum',wmin=.3,wstep=.2,nwin=40,level=.2,quietly=TRUE,approx=TRUE)
# plot(tmp$results[,7],tmp$results[,1],ylab='p-values',xlab='bandwidth')
# abline(v=1.3,lty='dashed')

# tmp <- rdwinselect(R,cbind(Plac[,1],X60),reps=wreps,stat='hotelling',wmin=.3,wstep=.2,nwin=40,level=.2,quietly=TRUE,approx=TRUE)
# plot(tmp$results[,7],tmp$results[,1],ylab='p-values',xlab='bandwidth')
# abline(v=2.7,lty='dashed')


###################################################################
## Table SA-9 and SA-10: Local Randomization Methods
###################################################################

# rreps <- 1000

# TableSA9_0 <- array(NA,dim=c(13,5))
# TableSA9_1 <- array(NA,dim=c(13,5))

# # p = 0

# tmp <- rdrandinf(Y,R,wl=-.9,wr=.9,reps=rreps,p=0)
# TableSA9_0[1,1] <- 0
# TableSA9_0[2,1] <- tmp$window[2]
# TableSA9_0[3,1] <- tmp$obs.stat
# TableSA9_0[4,1] <- tmp$p.value
# TableSA9_0[5,1] <- tmp$sumstats[2,1]
# TableSA9_0[6,1] <- tmp$sumstats[2,2]

# tmp <- rdrandinf(Y,R,wl=-1.1,wr=1.1,reps=rreps,p=0)
# TableSA9_0[1,2] <- 0
# TableSA9_0[2,2] <- tmp$window[2]
# TableSA9_0[3,2] <- tmp$obs.stat
# TableSA9_0[4,2] <- tmp$p.value
# TableSA9_0[5,2] <- tmp$sumstats[2,1]
# TableSA9_0[6,2] <- tmp$sumstats[2,2]

# tmp <- rdrandinf(Y,R,wl=-1.3,wr=1.3,reps=rreps,p=0)
# TableSA9_0[1,3] <- 0
# TableSA9_0[2,3] <- tmp$window[2]
# TableSA9_0[3,3] <- tmp$obs.stat
# TableSA9_0[4,3] <- tmp$p.value
# TableSA9_0[5,3] <- tmp$sumstats[2,1]
# TableSA9_0[6,3] <- tmp$sumstats[2,2]

# tmp <- rdrandinf(Y,R,wl=-1.5,wr=1.5,reps=rreps,p=0)
# TableSA9_0[1,4] <- 0
# TableSA9_0[2,4] <- tmp$window[2]
# TableSA9_0[3,4] <- tmp$obs.stat
# TableSA9_0[4,4] <- tmp$p.value
# TableSA9_0[5,4] <- tmp$sumstats[2,1]
# TableSA9_0[6,4] <- tmp$sumstats[2,2]

# tmp <- rdrandinf(Y,R,wl=-2.7,wr=2.7,reps=rreps,p=0)
# TableSA9_0[1,5] <- 0
# TableSA9_0[2,5] <- tmp$window[2]
# TableSA9_0[3,5] <- tmp$obs.stat
# TableSA9_0[4,5] <- tmp$p.value
# TableSA9_0[5,5] <- tmp$sumstats[2,1]
# TableSA9_0[6,5] <- tmp$sumstats[2,2]

# row <- 7
# for(col in 1:ncol(Plac)){
#   tmp <- rdrandinf(Plac[,col],R,wl=-.9,wr=.9,reps=rreps,p=0)
#   TableSA9_0[row,1] <- tmp$p.value
#   tmp <- rdrandinf(Plac[,col],R,wl=-1.1,wr=1.1,reps=rreps,p=0)
#   TableSA9_0[row,2] <- tmp$p.value
#   tmp <- rdrandinf(Plac[,col],R,wl=-1.3,wr=1.3,reps=rreps,p=0)
#   TableSA9_0[row,3] <- tmp$p.value
#   tmp <- rdrandinf(Plac[,col],R,wl=-1.5,wr=1.5,reps=rreps,p=0)
#   TableSA9_0[row,4] <- tmp$p.value
#   tmp <- rdrandinf(Plac[,col],R,wl=-2.7,wr=2.7,reps=rreps,p=0)
#   TableSA9_0[row,5] <- tmp$p.value
#   row <- row + 1
# }

# # p = 1

# tmp <- rdrandinf(Y,R,wl=-.9,wr=.9,reps=rreps,p=1)
# TableSA9_1[1,1] <- 1
# TableSA9_1[2,1] <- tmp$window[2]
# TableSA9_1[3,1] <- tmp$obs.stat
# TableSA9_1[4,1] <- tmp$p.value
# TableSA9_1[5,1] <- tmp$sumstats[2,1]
# TableSA9_1[6,1] <- tmp$sumstats[2,2]

# tmp <- rdrandinf(Y,R,wl=-1.1,wr=1.1,reps=rreps,p=1)
# TableSA9_1[1,2] <- 1
# TableSA9_1[2,2] <- tmp$window[2]
# TableSA9_1[3,2] <- tmp$obs.stat
# TableSA9_1[4,2] <- tmp$p.value
# TableSA9_1[5,2] <- tmp$sumstats[2,1]
# TableSA9_1[6,2] <- tmp$sumstats[2,2]

# tmp <- rdrandinf(Y,R,wl=-1.3,wr=1.3,reps=rreps,p=1)
# TableSA9_1[1,3] <- 1
# TableSA9_1[2,3] <- tmp$window[2]
# TableSA9_1[3,3] <- tmp$obs.stat
# TableSA9_1[4,3] <- tmp$p.value
# TableSA9_1[5,3] <- tmp$sumstats[2,1]
# TableSA9_1[6,3] <- tmp$sumstats[2,2]

# tmp <- rdrandinf(Y,R,wl=-1.5,wr=1.5,reps=rreps,p=1)
# TableSA9_1[1,4] <- 1
# TableSA9_1[2,4] <- tmp$window[2]
# TableSA9_1[3,4] <- tmp$obs.stat
# TableSA9_1[4,4] <- tmp$p.value
# TableSA9_1[5,4] <- tmp$sumstats[2,1]
# TableSA9_1[6,4] <- tmp$sumstats[2,2]

# tmp <- rdrandinf(Y,R,wl=-2.7,wr=2.7,reps=rreps,p=1)
# TableSA9_1[1,5] <- 1
# TableSA9_1[2,5] <- tmp$window[2]
# TableSA9_1[3,5] <- tmp$obs.stat
# TableSA9_1[4,5] <- tmp$p.value
# TableSA9_1[5,5] <- tmp$sumstats[2,1]
# TableSA9_1[6,5] <- tmp$sumstats[2,2]

# row <- 7
# for(col in 1:ncol(Plac)){
#   tmp <- rdrandinf(Plac[,col],R,wl=-.9,wr=.9,reps=rreps,p=1)
#   TableSA9_1[row,1] <- tmp$p.value
#   tmp <- rdrandinf(Plac[,col],R,wl=-1.1,wr=1.1,reps=rreps,p=1)
#   TableSA9_1[row,2] <- tmp$p.value
#   tmp <- rdrandinf(Plac[,col],R,wl=-1.3,wr=1.3,reps=rreps,p=1)
#   TableSA9_1[row,3] <- tmp$p.value
#   tmp <- rdrandinf(Plac[,col],R,wl=-1.5,wr=1.5,reps=rreps,p=1)
#   TableSA9_1[row,4] <- tmp$p.value
#   tmp <- rdrandinf(Plac[,col],R,wl=-2.7,wr=2.7,reps=rreps,p=1)
#   TableSA9_1[row,5] <- tmp$p.value
#   row <- row + 1
# }

# round(TableSA9_0,3)
# round(TableSA9_1,3)

# # Asymptotic

# TableSA10_0 <- array(NA,dim=c(13,5))
# TableSA10_1 <- array(NA,dim=c(13,5))

# # p = 0

# tmp <- rdrandinf(Y,R,wl=-.9,wr=.9,reps=rreps,p=0)
# TableSA10_0[1,1] <- 0
# TableSA10_0[2,1] <- tmp$window[2]
# TableSA10_0[3,1] <- tmp$obs.stat
# TableSA10_0[4,1] <- tmp$asy.pvalue
# TableSA10_0[5,1] <- tmp$sumstats[2,1]
# TableSA10_0[6,1] <- tmp$sumstats[2,2]

# tmp <- rdrandinf(Y,R,wl=-1.1,wr=1.1,reps=rreps,p=0)
# TableSA10_0[1,2] <- 0
# TableSA10_0[2,2] <- tmp$window[2]
# TableSA10_0[3,2] <- tmp$obs.stat
# TableSA10_0[4,2] <- tmp$asy.pvalue
# TableSA10_0[5,2] <- tmp$sumstats[2,1]
# TableSA10_0[6,2] <- tmp$sumstats[2,2]

# tmp <- rdrandinf(Y,R,wl=-1.3,wr=1.3,reps=rreps,p=0)
# TableSA10_0[1,3] <- 0
# TableSA10_0[2,3] <- tmp$window[2]
# TableSA10_0[3,3] <- tmp$obs.stat
# TableSA10_0[4,3] <- tmp$asy.pvalue
# TableSA10_0[5,3] <- tmp$sumstats[2,1]
# TableSA10_0[6,3] <- tmp$sumstats[2,2]

# tmp <- rdrandinf(Y,R,wl=-1.5,wr=1.5,reps=rreps,p=0)
# TableSA10_0[1,4] <- 0
# TableSA10_0[2,4] <- tmp$window[2]
# TableSA10_0[3,4] <- tmp$obs.stat
# TableSA10_0[4,4] <- tmp$asy.pvalue
# TableSA10_0[5,4] <- tmp$sumstats[2,1]
# TableSA10_0[6,4] <- tmp$sumstats[2,2]

# tmp <- rdrandinf(Y,R,wl=-2.7,wr=2.7,reps=rreps,p=0)
# TableSA10_0[1,5] <- 0
# TableSA10_0[2,5] <- tmp$window[2]
# TableSA10_0[3,5] <- tmp$obs.stat
# TableSA10_0[4,5] <- tmp$asy.pvalue
# TableSA10_0[5,5] <- tmp$sumstats[2,1]
# TableSA10_0[6,5] <- tmp$sumstats[2,2]

# row <- 7
# for(col in 1:ncol(Plac)){
#   tmp <- rdrandinf(Plac[,col],R,wl=-.9,wr=.9,reps=rreps,p=0)
#   TableSA10_0[row,1] <- tmp$asy.pvalue
#   tmp <- rdrandinf(Plac[,col],R,wl=-1.1,wr=1.1,reps=rreps,p=0)
#   TableSA10_0[row,2] <- tmp$asy.pvalue
#   tmp <- rdrandinf(Plac[,col],R,wl=-1.3,wr=1.3,reps=rreps,p=0)
#   TableSA10_0[row,3] <- tmp$asy.pvalue
#   tmp <- rdrandinf(Plac[,col],R,wl=-1.5,wr=1.5,reps=rreps,p=0)
#   TableSA10_0[row,4] <- tmp$asy.pvalue
#   tmp <- rdrandinf(Plac[,col],R,wl=-2.7,wr=2.7,reps=rreps,p=0)
#   TableSA10_0[row,5] <- tmp$asy.pvalue
#   row <- row + 1
# }

# # p = 1

# tmp <- rdrandinf(Y,R,wl=-.9,wr=.9,reps=rreps,p=1)
# TableSA10_1[1,1] <- 1
# TableSA10_1[2,1] <- tmp$window[2]
# TableSA10_1[3,1] <- tmp$obs.stat
# TableSA10_1[4,1] <- tmp$asy.pvalue
# TableSA10_1[5,1] <- tmp$sumstats[2,1]
# TableSA10_1[6,1] <- tmp$sumstats[2,2]

# tmp <- rdrandinf(Y,R,wl=-1.1,wr=1.1,reps=rreps,p=1)
# TableSA10_1[1,2] <- 1
# TableSA10_1[2,2] <- tmp$window[2]
# TableSA10_1[3,2] <- tmp$obs.stat
# TableSA10_1[4,2] <- tmp$asy.pvalue
# TableSA10_1[5,2] <- tmp$sumstats[2,1]
# TableSA10_1[6,2] <- tmp$sumstats[2,2]

# tmp <- rdrandinf(Y,R,wl=-1.3,wr=1.3,reps=rreps,p=1)
# TableSA10_1[1,3] <- 1
# TableSA10_1[2,3] <- tmp$window[2]
# TableSA10_1[3,3] <- tmp$obs.stat
# TableSA10_1[4,3] <- tmp$asy.pvalue
# TableSA10_1[5,3] <- tmp$sumstats[2,1]
# TableSA10_1[6,3] <- tmp$sumstats[2,2]

# tmp <- rdrandinf(Y,R,wl=-1.5,wr=1.5,reps=rreps,p=1)
# TableSA10_1[1,4] <- 1
# TableSA10_1[2,4] <- tmp$window[2]
# TableSA10_1[3,4] <- tmp$obs.stat
# TableSA10_1[4,4] <- tmp$asy.pvalue
# TableSA10_1[5,4] <- tmp$sumstats[2,1]
# TableSA10_1[6,4] <- tmp$sumstats[2,2]

# tmp <- rdrandinf(Y,R,wl=-2.7,wr=2.7,reps=rreps,p=1)
# TableSA10_1[1,5] <- 1
# TableSA10_1[2,5] <- tmp$window[2]
# TableSA10_1[3,5] <- tmp$obs.stat
# TableSA10_1[4,5] <- tmp$asy.pvalue
# TableSA10_1[5,5] <- tmp$sumstats[2,1]
# TableSA10_1[6,5] <- tmp$sumstats[2,2]

# row <- 7
# for(col in 1:ncol(Plac)){
#   tmp <- rdrandinf(Plac[,col],R,wl=-.9,wr=.9,reps=rreps,p=1)
#   TableSA10_1[row,1] <- tmp$asy.pvalue
#   tmp <- rdrandinf(Plac[,col],R,wl=-1.1,wr=1.1,reps=rreps,p=1)
#   TableSA10_1[row,2] <- tmp$asy.pvalue
#   tmp <- rdrandinf(Plac[,col],R,wl=-1.3,wr=1.3,reps=rreps,p=1)
#   TableSA10_1[row,3] <- tmp$asy.pvalue
#   tmp <- rdrandinf(Plac[,col],R,wl=-1.5,wr=1.5,reps=rreps,p=1)
#   TableSA10_1[row,4] <- tmp$asy.pvalue
#   tmp <- rdrandinf(Plac[,col],R,wl=-2.7,wr=2.7,reps=rreps,p=1)
#   TableSA10_1[row,5] <- tmp$asy.pvalue
#   row <- row + 1
# }

# round(TableSA10_0,3)
# round(TableSA10_1,3)


###################################################################
## Figure SA-6: Sensitivity of Linear Adjustment Model
###################################################################

# ii <- which(abs(R)<=1.1)
# il <- which(abs(R)<=1.1 & D==0)
# ir <- which(abs(R)<=1.1 & D==1)

# plot(R[ii],Y[ii])
# reg.l <- lm(Y[il]~R[il])
# clip(-1,0,0,80)
# abline(lm(Y[il]~R[il]))
# clip(0,1,0,80)
# abline(lm(Y[ir]~R[ir]))

# ii <- which(abs(R)<=1.3)
# il <- which(abs(R)<=1.3 & D==0)
# ir <- which(abs(R)<=1.3 & D==1)

# plot(R[ii],Y[ii])
# reg.l <- lm(Y[il]~R[il])
# clip(-1,0,0,80)
# abline(lm(Y[il]~R[il]))
# clip(0,1,0,80)
# abline(lm(Y[ir]~R[ir]))


###################################################################
## Table SA-11: Comparison of Inference Approaches
###################################################################

## NOTE: this table is generated at the end of the main paper replication do file.
## See Figure 4: Summary of Results.