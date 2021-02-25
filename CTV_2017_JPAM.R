################################################################################
## Comparing Inference Approaches for RD Designs:
## A Reexamination of the Effect of Head Start on Child Mortality
## Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
## Last update: 22-FEB-2021
################################################################################
## SOFTWARE WEBSITE: https://rdpackages.github.io/
################################################################################
## TO INSTALL/DOWNLOAD R PACKAGES/FUNCTIONS:
## RDROBUST: install.packages('rdrobust')
## RDDENSITY: install.packages('rddensity')
## RDLOCRAND: install.packages('rdlocrand')
################################################################################
## NOTE: if you are using RDROBUST version 2020 or newer, the option 
## masspoints="off" and stdvars="on" may be needed to replicate the results 
## in the paper.
## For example, line 279:
##    tmp = rdrobust(Y,R,p=0,q=1)
## should be replaced by:
##    tmp = rdrobust(Y,R,p=0,q=1,masspoints="off",stdvars="on")
################################################################################
## NOTE: if you are using RDDENSITY version 2020 or newer, the option 
## masspoints=FALSE may be needed to replicate the results in the paper.
## For example, line 118:
##    tmp = rddensity(Rraw,c=cutoff)
## should be replaced by:
##    tmp = rddensity(Rraw,c=cutoff,masspoints=FALSE)
################################################################################

rm(list = ls())
options(width=200)

library(rdlocrand)
library(rdrobust)
library(sandwich)
library(rddensity)

#############################################################
## Load and setup data
#############################################################

data <- read.csv("headstart.csv")

dim(data)
names(data)

# 1960 Census covariates

cutoff <- 59.1984

X60 <- cbind(data$census1960_pop,
            data$census1960_pctsch1417,
            data$census1960_pctsch534, 
            data$census1960_pctsch25plus,
            data$census1960_pop1417,
            data$census1960_pop534,
            data$census1960_pop25plus,
            data$census1960_pcturban,
            data$census1960_pctblack)

# 1990 Census covariates

X90 <- cbind(data$census1990_pop,
            data$census1990_pop1824,
            data$census1990_pop2534,
            data$census1990_pop3554,
            data$census1990_pop55plus,
            data$census1990_pcturban,
            data$census1990_pctblack,
            data$census1990_percapinc)

# Outcome, normalized running variable and treatment

Y <- data$mort_age59_related_postHS
Rraw <- data$povrate60
R <- data$povrate60 - cutoff

D <- as.numeric(R>=0)

# Placebo outcomes

Plac <- cbind(data$mort_age59_injury_postHS,data$mort_age59_related_preHS)


#############################################################
## Figure 1: Scatter and RD Plot, Head Start Data
#############################################################

ii <- which(Y<20 & !is.na(Y) & !is.na(Rraw))
rdplot(Y[ii],Rraw[ii],c=59.1984,nbins=3000)
rdplot(Y[ii],Rraw[ii],c=59.1984)

###################################################################
## Table 1: Binomial Tests
###################################################################

Table1 <- array(NA,dim=c(6,4))

row <- 1
for(w in seq(.3,1.3,by=.2)){
  ww <- which(abs(R)<=w & !is.na(R) & !is.na(Y))
  Table1 [row,1] <- w
  Table1 [row,2] <- sum(1-D[ww])
  Table1 [row,3] <- sum(D[ww])
  tmp <- binom.test(sum(D[ww]),length(D[ww]),p=.5)
  Table1[row,4] <- tmp$p.value
  row <- row + 1
}

round(Table1,3)

###################################################################
## Table 2: Nonparametric Density Continuity Tests
###################################################################

Table2 <- array(NA,dim=c(3,5))

tmp <- rddensity(Rraw,c=cutoff)
Table2[1,1] <- tmp$h[[1]]
Table2[1,2] <- tmp$h[[2]]
Table2[1,3] <- tmp$N[[4]]
Table2[1,4] <- tmp$N[[5]]
Table2[1,5] <- tmp$test[[4]]

tmp <- rddensity(Rraw,c=cutoff,bwselect='diff')
Table2[2,1] <- tmp$h[[1]]
Table2[2,2] <- tmp$h[[2]]
Table2[2,3] <- tmp$N[[4]]
Table2[2,4] <- tmp$N[[5]]
Table2[2,5] <- tmp$test[[4]]

tmp <- rddensity(Rraw,c=cutoff,fitselect='restricted')
Table2[3,1] <- tmp$h[[1]]
Table2[3,2] <- tmp$h[[2]]
Table2[3,3] <- tmp$N[[4]]
Table2[3,4] <- tmp$N[[5]]
Table2[3,5] <- tmp$test[[4]]

round(Table2,3)

#############################################################
## Table 3: Flexible Parametric RD Methods
#############################################################

Table3 <- array(NA,dim=c(10,4))

## Main outcomes

# Linear model

col <- 1
p <- 1

for(h in c(9,18)){
  ii <- which(abs(R)<=h)
  ir <- which(R<=h & R>=0)
  il <- which(R>=-h & R<0)
  
  Nr <- sum(R<=h & R>=0 & !is.na(R) & !is.na(Y))
  Nl <- sum(R>=-h & R<0 & !is.na(R) & !is.na(Y))
  
  P <- poly(R[ii],p,raw=TRUE)
  
  tmp <- lm(Y[ii]~D[ii]+P+I(P*D[ii]))
  
  beta <- tmp$coef['D[ii]']
  se.aux <- sqrt(diag(vcovHC(tmp,type='HC1')))
  se <- se.aux[2]
  tstat <- beta/se
  
  Table3[1,col] <- p
  Table3[2,col] <- h
  Table3[3,col] <- beta
  Table3[4,col] <- beta-qt(.975,Nr+Nl)*se
  Table3[5,col] <- beta+qt(.975,Nr+Nl)*se
  Table3[6,col] <- 2*(1-pt(abs(tstat),Nr+Nl))
  Table3[7,col] <- Nl
  Table3[8,col] <- Nr
  
  col <- col +1
}


# Quartic model

p <- 4

for(h in c(20,100)){
ii <- which(abs(R)<=h)
  ir <- which(R<=h & R>=0)
  il <- which(R>=-h & R<0)
  
  Nr <- sum(R<=h & R>=0 & !is.na(R) & !is.na(Y))
  Nl <- sum(R>=-h & R<0 & !is.na(R) & !is.na(Y))
  
  P <- poly(R[ii],p,raw=TRUE)
  
  tmp <- lm(Y[ii]~D[ii]+P+I(P*D[ii]))
  
  beta <- tmp$coef['D[ii]']
  se.aux <- sqrt(diag(vcovHC(tmp,type='HC1')))
  se <- se.aux[2]
  tstat <- beta/se
  
  Table3[1,col] <- p
  Table3[2,col] <- h
  Table3[3,col] <- beta
  Table3[4,col] <- beta-qt(.975,Nr+Nl)*se
  Table3[5,col] <- beta+qt(.975,Nr+Nl)*se
  Table3[6,col] <- 2*(1-pt(abs(tstat),Nr+Nl))
  Table3[7,col] <- Nl
  Table3[8,col] <- Nr
  
  col <- col + 1
}

## Placebo outcomes

# Linear model

col <- 1
p <- 1
for(h in c(9,18)){
  row <- 9
  ii <- which(abs(R)<=h)
  ir <- which(R<=h & R>=0)
  il <- which(R>=-h & R<0)
  P <- poly(R[ii],p,raw=TRUE)

  for(var in c(1,2)){
    Nr <- sum(R<=h & R>=0 & !is.na(R) & !is.na(Plac[,var]))
    Nl <- sum(R>=-h & R<0 & !is.na(R) & !is.na(Plac[,var]))
    tmp <- lm(Plac[ii,var]~D[ii]+P+I(P*D[ii]))
    beta <- tmp$coef['D[ii]']
    se.aux <- sqrt(diag(vcovHC(tmp,type='HC1')))
    se <- se.aux[2]
    tstat <- beta/se
    Table3[row,col] <- 2*(1-pt(abs(tstat),Nr+Nl))
    row <- row + 1
  }
  col <- col + 1
}


# Quartic model

p <- 4
for(h in c(20,100)){
  row <- 9
  ii <- which(abs(R)<=h)
  ir <- which(R<=h & R>=0)
  il <- which(R>=-h & R<0)
  P <- poly(R[ii],p,raw=TRUE)

  for(var in c(1,2)){
    Nr <- sum(R<=h & R>=0 & !is.na(R) & !is.na(Plac[,var]))
    Nl <- sum(R>=-h & R<0 & !is.na(R) & !is.na(Plac[,var]))
    tmp <- lm(Plac[ii,var]~D[ii]+P+I(P*D[ii]))
    beta <- tmp$coef['D[ii]']
    se.aux <- sqrt(diag(vcovHC(tmp,type='HC1')))
    se <- se.aux[2]
    tstat <- beta/se
    Table3[row,col] <- 2*(1-pt(abs(tstat),Nr+Nl))
    row <- row + 1
  }
  col <- col + 1
}

round(Table3,3)

#############################################################
## Table 4: Robust Nonparametric Local polynomial Methods
#############################################################

Table4 <- array(NA,dim=c(10,4))

# Local Constant Regression

tmp <- rdrobust(Y,R,p=0,q=1)
Table4[1,1] <- tmp$p
Table4[2,1] <- tmp$bws[1,1]
Table4[3,1] <- tmp$coef[1]
Table4[4,1] <- tmp$ci[3,1]
Table4[5,1] <- tmp$ci[3,2]
Table4[6,1] <- tmp$pv[3]
Table4[7,1] <- tmp$N_h[1]
Table4[8,1] <- tmp$N_h[2]


tmp <- rdrobust(Y,R,p=0,q=1,h=9)
Table4[1,2] <- tmp$p
Table4[2,2] <- tmp$bws[1,1]
Table4[3,2] <- tmp$coef[1]
Table4[4,2] <- tmp$ci[3,1]
Table4[5,2] <- tmp$ci[3,2]
Table4[6,2] <- tmp$pv[3]
Table4[7,2] <- tmp$N_h[1]
Table4[8,2] <- tmp$N_h[2]


tmp <- rdrobust(Plac[,1],R,p=0,q=1)
Table4[9,1] <- tmp$pv[3]
tmp <- rdrobust(Plac[,2],R,p=0,q=1)
Table4[10,1] <- tmp$pv[3]

tmp <- rdrobust(Plac[,1],R,p=0,q=1,h=9)
Table4[9,2] <- tmp$pv[3]
tmp <- rdrobust(Plac[,2],R,p=0,q=1,h=9)
Table4[10,2] <- tmp$pv[3]

# Local Linear Regression

tmp <- rdrobust(Y,R,p=1)
Table4[1,3] <- tmp$p
Table4[2,3] <- tmp$bws[1,1]
Table4[3,3] <- tmp$coef[1]
Table4[4,3] <- tmp$ci[3,1]
Table4[5,3] <- tmp$ci[3,2]
Table4[6,3] <- tmp$pv[3]
Table4[7,3] <- tmp$N_h[1]
Table4[8,3] <- tmp$N_h[2]

tmp <- rdrobust(Y,R,p=1,h=9)
Table4[1,4] <- tmp$p
Table4[2,4] <- tmp$bws[1,1]
Table4[3,4] <- tmp$coef[1]
Table4[4,4] <- tmp$ci[3,1]
Table4[5,4] <- tmp$ci[3,2]
Table4[6,4] <- tmp$pv[3]
Table4[7,4] <- tmp$N_h[1]
Table4[8,4] <- tmp$N_h[2]

tmp <- rdrobust(Plac[,1],R,p=1)
Table4[9,3] <- tmp$pv[3]
tmp <- rdrobust(Plac[,2],R,p=1)
Table4[10,3] <- tmp$pv[3]

tmp <- rdrobust(Plac[,1],R,p=1,h=9)
Table4[9,4] <- tmp$pv[3]
tmp <- rdrobust(Plac[,2],R,p=1,h=9)
Table4[10,4] <- tmp$pv[3]

round(Table4,3)

#############################################################
## Figure 2: Window Selection and Outcome of Interest
#############################################################

## Window selection

Xrdw <- cbind(data$mort_age59_related_preHS,X60)
tmp <- rdwinselect(R,Xrdw,reps=1000,statistic="ksmirnov",wmin=.3,wstep=.2,level=.2)

## P-values plot

tmp <- rdwinselect(R,Xrdw,reps=1000,statistic="ksmirnov",wmin=.3,wstep=.2,level=.2,nwin=40,plot=TRUE,quietly=TRUE)

## Scatter plot with means

w <- 1.1
ir <- which(abs(R)<=w & D==1 & !is.na(Y) & !is.na(R))
il <- which(abs(R)<=w & D==0 & !is.na(Y) & !is.na(R))
ii <- which(abs(R)<=w & !is.na(Y) & !is.na(R))

ml <- mean(Y[il])
mr <- mean(Y[ir])

plot(R[ii],Y[ii])
segments(-w,ml,0,ml,lty='dashed')
segments(0,mr,w,mr,lty='dashed')

#############################################################
## Table 5: Local-Randomization Methods
#############################################################

## Inference results

w <- 1.1
reps <- 1000

Table5 <- array(NA,dim=c(8,6))

# Outcome

tmp <- rdrandinf(Y,R,wl=-w,wr=w,reps=reps,quietly=TRUE)
Table5[1:6,1] <- c(0,tmp$window[2],tmp$obs.stat,tmp$p.value,tmp$sumstats[2,1],tmp$sumstats[2,2])

tmp <- rdrandinf(Y,R,wl=-3.235,wr=3.235,reps=reps,quietly=TRUE)
Table5[1:6,2] <- c(0,tmp$window[2],tmp$obs.stat,tmp$p.value,tmp$sumstats[2,1],tmp$sumstats[2,2])

tmp <- rdrandinf(Y,R,wl=-9,wr=9,reps=reps,quietly=TRUE)
Table5[1:6,3] <- c(0,tmp$window[2],tmp$obs.stat,tmp$p.value,tmp$sumstats[2,1],tmp$sumstats[2,2])

tmp <- rdrandinf(Y,R,wl=-w,wr=w,reps=reps,p=1,quietly=TRUE)
Table5[1:6,4] <- c(1,tmp$window[2],tmp$obs.stat,tmp$p.value,tmp$sumstats[2,1],tmp$sumstats[2,2])

tmp <- rdrandinf(Y,R,wl=-3.235,wr=3.235,reps=reps,p=1,quietly=TRUE)
Table5[1:6,5] <- c(1,tmp$window[2],tmp$obs.stat,tmp$p.value,tmp$sumstats[2,1],tmp$sumstats[2,2])

tmp <- rdrandinf(Y,R,wl=-9,wr=9,reps=reps,p=1,quietly=TRUE)
Table5[1:6,6] <- c(1,tmp$window[2],tmp$obs.stat,tmp$p.value,tmp$sumstats[2,1],tmp$sumstats[2,2])

# Placebo outcomes

for (i in 1:ncol(Plac)){
  tmp <- rdrandinf(Plac[,i],R,wl=-w,wr=w,reps=reps,quietly=TRUE)
  Table5[6+i,1] <- tmp$p.value
  tmp <- rdrandinf(Plac[,i],R,wl=-3.235,wr=3.235,reps=reps,quietly=TRUE)
  Table5[6+i,2] <- tmp$p.value
  tmp <- rdrandinf(Plac[,i],R,wl=-9,wr=9,reps=reps,quietly=TRUE)
  Table5[6+i,3] <- tmp$p.value
  tmp <- rdrandinf(Plac[,i],R,wl=-w,wr=w,reps=reps,p=1,quietly=TRUE)
  Table5[6+i,4] <- tmp$p.value
  tmp <- rdrandinf(Plac[,i],R,wl=-3.235,wr=3.235,reps=reps,p=1,quietly=TRUE)
  Table5[6+i,5] <- tmp$p.value
  tmp <- rdrandinf(Plac[,i],R,wl=-9,wr=9,reps=reps,p=1,quietly=TRUE)
  Table5[6+i,6] <- tmp$p.value
}

round(Table5,3)

#############################################################
## Figure 3: Sensitivity to Window Length
#############################################################

tmp <- rdsensitivity(Y,R,wlist=seq(.3,10.1,by=.2),tlist=seq(-10,6.5,by=.25))
tmp <- rdsensitivity(Y,R,wlist=seq(.3,10.1,by=.2),tlist=seq(-10,6.5,by=.25),p=1)

#############################################################
## Table 6: Local Randomization methods -- CI and interf
#############################################################

Table6 <- array(NA,dim=c(10,2))

reps <- 5000

tmp <- rdrandinf(Y,R,wl=-w,wr=w,reps=reps,ci=c(.05,seq(-5,0,by=.025)),quietly=TRUE)
Table6[c(1:6,9,10),1] <- c(0,tmp$window[2],tmp$obs.stat,tmp$p.value,tmp$ci[1],tmp$ci[2],tmp$sumstats[2,1],tmp$sumstats[2,2])

tmp <- rdrandinf(Y,R,wl=-w,wr=w,reps=reps,interfci=.05,quietly=TRUE)
Table6[c(7,8),1] <- c(tmp$interf.ci[1],tmp$interf.ci[2])

tmp <- rdrandinf(Y,R,wl=-w,wr=w,reps=reps,ci=c(.05,seq(-5,0,by=.025)),p=1,quietly=TRUE)
Table6[c(1:6,9,10),2] <- c(1,tmp$window[2],tmp$obs.stat,tmp$p.value,tmp$ci[1],tmp$ci[2],tmp$sumstats[2,1],tmp$sumstats[2,2])

tmp <- rdrandinf(Y,R,wl=-w,wr=w,reps=reps,interfci=.05,p=1,quietly=TRUE)
Table6[c(7,8),2] <- c(tmp$interf.ci[1],tmp$interf.ci[2])

round(Table6,3)

#############################################################
## Table 7: Rosenbaum Bounds
#############################################################

Table7 <- array(NA,dim=c(7,9))

tmp <- rdrbounds(Y,R,expgamma=c(1.1,1.2,1.3,1.4),wlist=seq(.3,1.5,by=.2),statistic='ttest',bound='upper',fmpval=TRUE,reps=5000)
Table7[1,3:9] <- seq(.3,1.5,by=.2)
Table7[2:3,3:9] <- tmp$p.values
Table7[4:7,3:9] <- tmp$upper.bound
Table7[4:7,1] <- log(c(1.1,1.2,1.3,1.4))
Table7[4:7,2] <- c(1.1,1.2,1.3,1.4)

round(Table7,3)

#############################################################
## Figure 4: Summary of Results
#############################################################

h <- c(1.1,1.3,3.235,6.811,9,18,25)
point.est <- c(Table5[3,1],Table5[3,4],Table4[3,3:4],Table3[3,1:2],Table3[3,4])
ci.lb <- c(Table6[5,1],Table6[5,2],Table4[4,3:4],Table3[4,1:2],Table3[4,4])
ci.rb <- c(Table6[6,1],Table6[6,2],Table4[5,3:4],Table3[5,1:2],Table3[5,4])

df <- data.frame(cbind(h,point.est))

library(ggplot2)

figure4 <- ggplot(df, aes(x=h,y=point.est)) + geom_point(size=4) + geom_errorbar(aes(ymax=ci.rb,ymin=ci.lb)) + ylim(-6,2)
figure4 <- figure4 + geom_hline(aes(yintercept=0),linetype='dashed') + geom_vline(aes(xintercept=2.2285)) + geom_vline(aes(xintercept=7.699)) + geom_vline(aes(xintercept=21.5))
figure4

#############################################################
## Additional Empirical Analysis
#############################################################

## Note: these results are not reported in the paper.

## Robust Nonparametric Method with Covariates

summary(rdrobust(Y,R,covs=X60))

##  Robust Nonparametric Method: Different Bandwdiths at Each Side

summary(rdrobust(Y,R,bwselect='msetwo'))
