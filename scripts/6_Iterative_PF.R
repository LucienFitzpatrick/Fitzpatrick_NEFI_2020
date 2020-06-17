#----------------------------------------------------------------------------------------------------------------------------------#
# Script by : NEFI 2020 ummer group
# Project: Living Collections Phenology Forecasting
# Purpose: To use arb weather data and phenology monitoring data to create a predicitve model of fall color intensity
#          This script serves as the framework for iterative forecasting
# Inputs: dat.npn dataframe containing 2019 npn observation data
#         A matrix output of a model containing variables of interest
# Outputs: Iterative hindcast of 2018 based on a model trained with 2019 data
# Notes: 
#-----------------------------------------------------------------------------------------------------------------------------------#

#####################################
# Particle Filter
library(ecoforecastR)

#This section is for defining a range for a 2018 hindcast
dat.2018 <- dat.2018[dat.2018$day_of_year > 213, ]
dat.2018$day_of_year <- as.numeric(as.character(dat.2018$day_of_year))
dat.2018$color.full <- as.numeric(as.character(dat.2018$color.full))

data_2018 <- list(y = dat.2018$color.full, n = length(dat.2018$color.full), time = dat.2018$day_of_year-213, nt = 365-213, 
                  a_add=0.1, r_add=0.001, x_ic = -100, tau_ic = 1000)

#Number of monte Carlo iterations that will be used
Nmc = 10000

#Set the Initial conditions using the parameters from your model
#IC <- rlnorm(Nmc, data$x_ic,(1/sqrt(data$tau_ic)))

#This if for 2018 hindcast testing
IC <- rlnorm(Nmc, data_2018$x_ic,(1/sqrt(data_2018$tau_ic)))

NT = data$nt
forecastx <- function(IC,Q=0,n=Nmc,betaLOD,LOD){
  x <- z <- matrix(NA,n,NT)  ## storage
  Xprev <- IC           ## initialize
  for(t in 1:NT){
    mu <- Xprev  + betaLOD*LOD[t] 
    z[,t] = rnorm(n,mu,Q)
    x[,t] <- pmin(0.999,pmax(Xprev,z[,t]))
    Xprev <- x[,t]                                  ## update IC
  }
  return(x)
}

param.mean <- apply(out_LOD,2,mean)

x.det_LOD <- forecastx(IC=mean(IC),
                       Q=mean((1/sqrt(out_LOD[,"tau_add"]))),  ## process error off
                       betaLOD=param.mean["betaLOD"],
                       LOD=dayLengths,
                       n=Nmc)

x.i_LOD <- forecastx(IC=IC,
                     Q=mean((1/sqrt(out_LOD[,"tau_add"]))),  ## process error off
                     betaLOD=param.mean["betaLOD"],
                     LOD=dayLengths,
                     n=Nmc)

s = seq(1, nrow(out_LOD), length=Nmc)

x.ip_LOD <- forecastx(IC=IC,
                      Q=1/sqrt(out_LOD[s,"tau_add"]),  ## process error off
                      betaLOD=out_LOD[s,"betaLOD"],
                      LOD=dayLengths,
                      n=Nmc)


det_ci_LOD <- apply(x.det_LOD,2,quantile,c(0.025,0.5,0.975)) 
i_ci_LOD <- apply(x.i_LOD,2,quantile,c(0.025,0.5,0.975)) 
ip_ci_LOD <- apply(x.ip_LOD,2,quantile,c(0.025,0.5,0.975)) 

## calculate the cumulative likelihoods
## to be used as PF weights
##Two matrices are used because the second is used for calculation, and then the second is filled for the correect time period
Npnlike = matrix(NA,nrow = Nmc, ncol = nrow(dat.2018))
Npnclike = matrix(0,nrow = Nmc, ncol = 365-213)
# i=1
for(i in 1:Nmc){
  Npnlike[i,] = dbinom(dat.2018$color.full,1,x.ip_LOD[i,dat.2018$day_of_year-213],log=TRUE)  ## calculate log likelihoods
  Npnlike[i,is.na(Npnlike[i,])] = 0       ## missing data as weight 1; log(1)=0
  tmp = tapply(Npnlike[i,],dat.2018$day_of_year-213,sum)
  tmp2 = rep(0, 365-213)
  tmp2[as.numeric(names(tmp))] = tmp
  Npnclike[i,] = exp(cumsum(tmp2))
}
hist(Npnclike[,ncol(Npnclike)],main="Final Ensemble Weights")

## Non-resampling Particle Filter
## calculation of CI
nobs = ncol(Npnclike)                     ## number of observations
Npnpf = matrix(NA,3,nobs)
wbar = apply(Npnclike,2,mean)             ## mean weight at each time point
for(i in 1:nobs){
  Npnpf[,i] = wtd.quantile(x.ip_LOD[,i],Npnclike[,i]/wbar[i],c(0.025,0.5,0.975))  ## calculate weighted median and CI
}
# Loess smoothing for 2018 data 
dat.2018_order = dat.2018[order(dat.2018$day_of_year),]
data_2018_loess_10 <- loess(as.numeric(as.character(color.full)) ~ as.numeric(as.character(day_of_year)), data=dat.2018_order, span=0.10) # 10% smoothing span
data_2018_loess_10 <- predict(data_2018_loess_10) 
data_2018_loess_30 <- loess(as.numeric(as.character(color.full)) ~ as.numeric(as.character(day_of_year)), data=dat.2018_order, span=0.30) # 10% smoothing span
data_2018_loess_30 <- predict(data_2018_loess_30) 


#Initial plotting of data and uncertainty partitioning BEFORE iteration is added
plot(dat.2018$day_of_year, dat.2018$color.full ,pch="+",cex=0.5, xlab = "Day of Year", ylab = "Fall Color", main = "2018 Non-resampling Particle Filter")
ecoforecastR::ciEnvelope(time,ip_ci_LOD[1,],ip_ci_LOD[3,],col=col.alpha("purple",0.5))
ecoforecastR::ciEnvelope(time,i_ci_LOD[1,],i_ci_LOD[3,],col=col.alpha("green",0.5))
ecoforecastR::ciEnvelope(time,det_ci_LOD[1,],det_ci_LOD[3,],col=col.alpha("blue",0.5))
ecoforecastR::ciEnvelope(time,Npnpf[1,],Npnpf[3,],col=col.alpha("red",0.5))
lines(data_2018_loess_10, x=dat.2018_order$day_of_year, col="green", lwd = 2)
lines(data_2018_loess_30, x=dat.2018_order$day_of_year, col="blue", lwd = 2)


s = seq(234, 355, by = 20) # date that we're running the forecast on 
# set weights after the date to 0

iter_pf = list()
i=1
for(i in seq_along(s)) {
  # what row index corresponds to date for forecasting up to
  iter_pf[[i]] = matrix(NA,3,nobs)
  wbar = apply(Npnclike,2,mean)             ## mean weight at each time point
  tmp_clike = Npnclike
  tmp_clike[,(s[i]:365)-213] <- tmp_clike[,s[i]-213]
  for(j in 1:nobs){
    iter_pf[[i]][,j] = wtd.quantile(x.ip_LOD[,j],tmp_clike[,j]/wbar[j],c(0.025,0.5,0.975))  ## calculate weighted median and CI
  }
}

plot(dat.2018$day_of_year, dat.2018$color.full ,pch="+",cex=0.5, xlab = "Day of Year", ylab = "Fall Color", main = "Iterative 2018 Non-resampling Particle Filter")
ecoforecastR::ciEnvelope(time,ip_ci_LOD[1,],ip_ci_LOD[3,],col=col.alpha("purple",0.5))
for(i in seq_along(s)) {
  ecoforecastR::ciEnvelope(time[(s[i]:365)-213],iter_pf[[i]][1,(s[i]:365)-213],iter_pf[[i]][3,(s[i]:365)-213],col=col.alpha(i,0.5))
}
# ecoforecastR::ciEnvelope(time,Npnpf[1,],Npnpf[3,],col=col.alpha("red",0.5))
lines(data_2018_loess_10, x=dat.2018_order$day_of_year, col="green", lwd = 2)
lines(data_2018_loess_30, x=dat.2018_order$day_of_year, col="blue", lwd = 2)
