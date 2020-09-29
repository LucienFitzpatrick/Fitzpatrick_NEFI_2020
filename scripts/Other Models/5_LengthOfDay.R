#----------------------------------------------------------------------------------------------------------------------------------#
# Script by : 
# Project: Living Collections Phenology Forecasting Autumn
# Purpose: To use arb weather data and phenology monitoring data to create a predicitve model of fall color intensity
#          This script is the state space model using day length as a covariate
# Inputs: Data frame of clean NPN observations
#         Data frame of daymet data
# Outputs: A model output R data object
# Notes: 
#-----------------------------------------------------------------------------------------------------------------------------------#

library(rjags)

path.doc <- ("../data_processed/fall/")

dat.npn <- read.csv(file.path(path.doc, "Fall_Phenology_data.csv"))


#creating 2018 frame for hindcasting
dat.2018 <- dat.npn[dat.npn$year == 2018, ]

dat.npn$color.clean <- as.numeric(as.character(dat.npn$color.clean))

time <- 214:365


# Model 
binom_LOD = "
model{

#### Data Model
for(i in 1:n){
y[i] ~ dbern(x[time[i]])
}

#### Process Model
for(t in 2:nt){
z[t]~dnorm(mu[t],tau_add)
mu[t] <- x[t-1]  + betaLOD*LOD[t] # consider: betaX*x[t-1] + betaIntercept 
x[t] <- min(0.999,max(0.0001,z[t]))
#x[t]~ dnorm(mu[t],tau_add)
}

#### Priors
x[1] ~ dlnorm(x_ic,tau_ic)
# tau_obs ~ dgamma(a_obs,r_obs)
tau_add ~ dgamma(a_add,r_add)
betaLOD ~ dnorm(0, 1)
}
"

data <- list(y = dat.npn$color.clean, n = length(dat.npn$color.clean), time = dat.npn$day_of_year-212, LOD = dat.npn$LOD, nt = 365-213, 
             a_add=1, r_add=1, x_ic =-100, tau_ic = 1000)




j.model_LOD   <- jags.model (file = textConnection(binom_LOD),
                         data = data,
                         n.chains = 3)

jags.out_LOD   <- coda.samples (model = j.model_LOD ,
                            variable.names = c("x","tau_add", "betaLOD"),
                            n.iter = 20000)

# GBR_LOD <- gelman.plot(jags.out_LOD)
gelman.diag(jags.out_LOD) # tau_add high
gelman.plot(jags.out_LOD[,"betaLOD"])

# rank order should be consistent, plot one that's the slowest to converge  

burnin_LOD = 10000                                ## determine convergence
jags.burn_LOD <- window(jags.out_LOD,start=burnin_LOD)  ## remove burn-in
#plot(jags.burn) 
summary(jags.burn_LOD)

out_LOD <- as.matrix(jags.burn_LOD)
path.hub <- "C:/Users/lucie/Documents/GitHub/NEFI/"
saveRDS(out_LOD, file.path(path.hub, "model_output/LengthOfDay_Output.RDS"))

#---------------------------------------------------------------------------------#
#This should be all you need in order to run the iterative particle filter script
#---------------------------------------------------------------------------------#

x.cols_LOD <- grep("^x",colnames(out_LOD)) ## grab all columns that start with the letter x
ci_LOD <- apply(out_LOD[,x.cols_LOD],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

time <- 214:365
data_all = dat.npn[,c("day_of_year","color.clean")]
data_all = data_all[order(data_all$day_of_year),]
data_all$color.clean = as.numeric(as.character(data_all$color.clean))
data_all$day_of_year = as.numeric(as.character(data_all$day_of_year))

data_all_loess_10 <- loess(as.numeric(as.character(color.clean)) ~ as.numeric(as.character(day_of_year)), data=data_all, span=0.10) # 10% smoothing span
data_all_loess_10 <- predict(data_all_loess_10) 
data_all_loess_30 <- loess(as.numeric(as.character(color.clean)) ~ as.numeric(as.character(day_of_year)), data=data_all, span=0.30) # 10% smoothing span
data_all_loess_30 <- predict(data_all_loess_30) 
plot(time,ci_LOD[2,],ylim=c(0,1),ylab="Fall color")
lines(data_all_loess_10, x=data_all$day_of_year, col="red", lwd = 2)
lines(data_all_loess_30, x=data_all$day_of_year, col="blue", lwd = 2)


plot(data_all$day_of_year, data_all$color.clean, type="p", main="Loess Smoothing 2019 Data", xlab="Day of Year", ylab="Fall Color")

## adjust x-axis label to be monthly if zoomed
ecoforecastR::ciEnvelope(time,ci_LOD[1,],ci_LOD[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(dat.npn$day_of_year, dat.npn$color.clean ,pch="+",cex=0.5)
lines(data_all_loess_10, x=data_all$day_of_year, col="red", lwd = 2)
lines(data_all_loess_30, x=data_all$day_of_year, col="blue", lwd = 2)

#dic_LOD = dic.samples(j.model_LOD, 5000)
#saveRDS(dic_LOD, "model_output/LengthOfDay_DIC.RDS")



NT = data$nt
forecastx <- function(IC,Q=0,n=Nmc,betaLOD,LOD){
  x <- z <- matrix(NA,n,NT)  ## storage
  Xprev <- IC           ## initialize
  for(t in 1:NT){
    mu <- Xprev  + betaLOD*LOD[t] 
    z[,t] = rnorm(n,mu,Q)
    x[,t] <- pmin(0.999,pmax(.0001,z[,t]))
    Xprev <- x[,t]                                  ## update IC
  }
  return(x)
}




Nmc = 10000
IC <- rlnorm(Nmc, data$x_ic,(1/sqrt(data$tau_ic)))
# IC <- as.matrix(out$predict)
# IC <- as.matrix(out_LOD[,3:154])


## calculate mean of all inputs
# ppt.mean <- matrix(apply(ppt_ensemble,2,mean),1,NT) ## driver
## parameters
# params <- as.matrix(jags.burn_LOD$params)
param.mean <- apply(out_LOD,2,mean)

x.det_LOD <- forecastx(IC=mean(IC),
                   Q=mean((1/sqrt(out_LOD[,"tau_add"]))),  ## process error off
                   betaLOD=param.mean["betaLOD"],
                   LOD=dat.npn$LOD,
                   n=Nmc)

x.i_LOD <- forecastx(IC=IC,
                 Q=mean((1/sqrt(out_LOD[,"tau_add"]))),  ## process error off
                 betaLOD=param.mean["betaLOD"],
                 LOD=dat.npn$LOD,
                 n=Nmc)

s = seq(1, nrow(out_LOD), length=Nmc)

x.ip_LOD <- forecastx(IC=IC,
                  Q=1/sqrt(out_LOD[s,"tau_add"]),  ## process error off
                  betaLOD=out_LOD[s,"betaLOD"],
                  LOD=dat.npn$LOD,
                  n=Nmc)

det_ci_LOD <- apply(x.det_LOD,2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
i_ci_LOD <- apply(x.i_LOD,2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
ip_ci_LOD <- apply(x.ip_LOD,2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale


plot(dat.npn$day_of_year, dat.npn$color.clean ,pch="+",cex=0.5, xlab = "Day of Year", ylab = "Fall Color", main = "2019")
library(ecoforecastR)
ecoforecastR::ciEnvelope(time,ip_ci_LOD[1,],ip_ci_LOD[3,],col=col.alpha("purple",0.5))
ecoforecastR::ciEnvelope(time,i_ci_LOD[1,],i_ci_LOD[3,],col=col.alpha("green",0.5))
ecoforecastR::ciEnvelope(time,det_ci_LOD[1,],det_ci_LOD[3,],col=col.alpha("blue",0.5))



# 2018 forecast
dat.2018
dat.2018 <- dat.2018[dat.2018$day_of_year > 213, ]
dat.2018$day_of_year <- as.numeric(as.character(dat.2018$day_of_year))
dat.2018$color.clean <- as.numeric(as.character(dat.2018$color.clean))

data_2018 <- list(y = dat.2018$color.clean, n = length(dat.2018$color.clean), time = dat.2018$day_of_year-213, nt = 365-213, 
             a_add=0.1, r_add=0.001, x_ic = -100, tau_ic = 1000)

Nmc = 10000
IC <- rlnorm(Nmc, data_2018$x_ic,(1/sqrt(data_2018$tau_ic)))

param.mean <- apply(out_LOD,2,mean)

x.det_LOD <- forecastx(IC=mean(IC),
                       Q=mean((1/sqrt(out_LOD[,"tau_add"]))),  ## process error off
                       betaLOD=param.mean["betaLOD"],
                       LOD=dat.npn$LOD,
                       n=Nmc)

x.i_LOD <- forecastx(IC=IC,
                     Q=mean((1/sqrt(out_LOD[,"tau_add"]))),  ## process error off
                     betaLOD=param.mean["betaLOD"],
                     LOD=dat.npn$LOD,
                     n=Nmc)

s = seq(1, nrow(out_LOD), length=Nmc)

x.ip_LOD <- forecastx(IC=IC,
                      Q=1/sqrt(out_LOD[s,"tau_add"]),  ## process error off
                      betaLOD=out_LOD[s,"betaLOD"],
                      LOD=dat.npn$LOD,
                      n=Nmc)

det_ci_LOD <- apply(x.det_LOD,2,quantile,c(0.025,0.5,0.975)) 
i_ci_LOD <- apply(x.i_LOD,2,quantile,c(0.025,0.5,0.975)) 
ip_ci_LOD <- apply(x.ip_LOD,2,quantile,c(0.025,0.5,0.975)) 


plot(dat.2018$day_of_year, dat.2018$color.clean ,pch="+",cex=0.5, xlab = "Day of Year", ylab = "Fall Color", main = "2018 Forecast")
ecoforecastR::ciEnvelope(time,ip_ci_LOD[1,],ip_ci_LOD[3,],col=col.alpha("purple",0.5))
ecoforecastR::ciEnvelope(time,i_ci_LOD[1,],i_ci_LOD[3,],col=col.alpha("green",0.5))
ecoforecastR::ciEnvelope(time,det_ci_LOD[1,],det_ci_LOD[3,],col=col.alpha("blue",0.5))
lines(data_2018_loess_10, x=dat.2018_order$day_of_year, col="red", lwd = 2)
lines(data_2018_loess_30, x=dat.2018_order$day_of_year, col="blue", lwd = 2)
