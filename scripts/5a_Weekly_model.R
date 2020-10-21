#----------------------------------------------------------------------------------------------------------------------------------#
# Script by : Lucien Fitzpatrick
# Project: Living Collections Phenology Forecasting
# Purpose: To use arb weather data and phenology monitoring data to create a predicitve model of fall color intensity
#          This script is the first, most simplified, version of the model that we have
# Inputs: Clean NPN observation dataframe from 3_NPN_clean.R script
#         Daymet_clean_data from 2_Daymet_download.R
# Outputs: Dataframe that has the weather calculations for all days of interest
# Notes: 
#-----------------------------------------------------------------------------------------------------------------------------------#
library(rjags)
library(coda)
library(ecoforecastR)


path.doc <- ("../data_processed/fall/")
path.fig <- ("../data_processed/fall/figures")

dat.npn <- read.csv(file.path(path.doc, "Fall_Phenology_data.csv"))

dat.npn$color.clean <- as.numeric(as.character(dat.npn$color.clean))


RandomWalk_binom = "
model{
  
  #### Data Model
  for(i in 1:n){
    y[i] ~ dbern(x[time[i]])
  }
  
  #### Process Model
  for(t in 2:nt){
    z[t] ~ dnorm(x[t-1],tau_add)
    x[t] <- min(0.999,max(0.0001,z[t]))
  }
  
  #### Priors
  x[1] ~ dlnorm(x_ic,tau_ic)
  tau_add ~ dgamma(a_add,r_add)
}
"


dat.npn <- dat.npn[order(dat.npn$week),]

#Removing weeks with less than .75 proportion of trees observed
dat.new <- dat.npn[dat.npn$week < 48,] 

data.new <- list(y = dat.new$color.clean, n = length(dat.new$color.clean), time = dat.new$week-(min(dat.new$week)-1), nt = length(unique(dat.new$week)), 
                 a_add=100, r_add=1, x_ic = -100 , tau_ic = 1000)



j.model   <- jags.model (file = textConnection(RandomWalk_binom),
                         data = data.new,
                         n.chains = 3)


jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add"),
                            n.iter = 40000)


gelman.diag(jags.out)

#plot(jags.out)

#GBR <- gelman.plot(jags.out)

burnin = 20000                                ## determine convergence
jags.burn <- window(jags.out,start=burnin)  ## remove burn-in
#plot(jags.burn)  
summary(jags.burn)## check diagnostics post burn-in

out <- as.matrix(jags.out)
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975))

time <- 30:47

data_all = dat.new[,c("week","color.clean")]
data_all = data_all[order(data_all$week),]
data_all$color.clean = as.numeric(as.character(data_all$color.clean))
data_all$week = as.numeric(as.character(data_all$week))

data_all_loess_10 <- loess(as.numeric(as.character(color.clean)) ~ as.numeric(as.character(week)), data=data_all, span=0.10) # 10% smoothing span
data_all_loess_10 <- predict(data_all_loess_10) 
data_all_loess_30 <- loess(as.numeric(as.character(color.clean)) ~ as.numeric(as.character(week)), data=data_all, span=0.30) # 30% smoothing span
data_all_loess_30 <- predict(data_all_loess_30) 


png(filename= file.path(path.fig, paste0("Weekly_Oak_Collection_Model_CI.png")))
plot(time,ci[2,],ylim=c(0,1),main ="Weekly_Oak_Collection_Model_CI.png", ylab="Fall color")
ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(dat.new$week, dat.new$color.clean ,pch="+",cex=0.5)
lines(data_all_loess_10, x=data_all$week, col="red", lwd = 2)
lines(data_all_loess_30, x=data_all$week, col="blue", lwd = 2)
dev.off()


#This section is for defining a range for a 2018 hindcast
dat.2018 <- dat.npn[dat.npn$year == 2018, ]
dat.2018 <- dat.2018[dat.2018$week > 29 & dat.2018$week< 48, ]
dat.2018$week <- as.numeric(as.character(dat.2018$week))
dat.2018$color.clean <- as.numeric(as.character(dat.2018$color.clean))

data_2018 <- list(y = dat.2018$color.clean, n = length(dat.2018$color.clean), time = dat.2018$week-(min(dat.2018$week)-1), nt = length(unique(dat.2018$week)), 
                  a_add=100, r_add=1, x_ic = -100, tau_ic = 1000)


#------------------------------------------------------------------------------------------#
#Here is where the 2018 forecast begins and things start getting tricky
#------------------------------------------------------------------------------------------#


#-------------------------------------#
#First is the forecast function that determines our uncertainty at each time point
#--------------------------------------#


#Number of monte Carlo iterations that will be used
Nmc = 20000

#Set the Initial conditions using the parameters from your model
#IC <- rlnorm(Nmc, data$x_ic,(1/sqrt(data$tau_ic)))

#This if for 2018 hindcast testing
IC <- rlnorm(Nmc, data_2018$x_ic,(1/sqrt(data_2018$tau_ic)))

#### Forecast function 
##` @param IC    Initial Conditions
##` @param Q     Process error (default = 0 for deterministic runs)
##` @param n     Size of Monte Carlo ensemble
NT = data_2018$nt
forecastx <- function(IC,Q=0,n=Nmc){
  x <- z <- matrix(NA,n,NT)  ## storage
  Xprev <- IC           ## initialize
  for(t in 1:NT){
    mu <- Xprev 
    z[,t] = rnorm(n,mu,Q)
    x[,t] <- pmin(0.999,pmax(0.0001,z[,t]))
    Xprev <- x[,t]                                  ## update IC
  }
  return(x)
}

param.mean <- apply(out,2,mean)

x.det_LOD <- forecastx(IC=mean(IC),
                       Q=mean((1/sqrt(out[,"tau_add"]))),  ## process error off
                       n=Nmc)

x.i_LOD <- forecastx(IC=IC,
                     Q=mean((1/sqrt(out[,"tau_add"]))),  ## process error off
                     n=Nmc)

s = seq(1, nrow(out), length=Nmc)

x.ip_LOD <- forecastx(IC=IC,
                      Q=1/sqrt(out[s,"tau_add"]),  ## process error off
                      n=Nmc)


det_ci_LOD <- apply(x.det_LOD,2,quantile,c(0.025,0.5,0.975)) 
i_ci_LOD <- apply(x.i_LOD,2,quantile,c(0.025,0.5,0.975)) 
ip_ci_LOD <- apply(x.ip_LOD,2,quantile,c(0.025,0.5,0.975)) 


#-----------------------------------------------------------------------------#
#Here is where the errors actually begin (in the rolling model)
#-----------------------------------------------------------------------------#

## calculate the cumulative likelihoods
## to be used as PF weights
#Two tmp are used because the first is used for calculation, and then the second is filled for the correct time period
Npnlike = matrix(NA,nrow = Nmc, ncol = nrow(dat.2018))#log liklihood
Npnclike = matrix(0,nrow = Nmc, ncol = length(unique(dat.2018$week)))#cumulative log liklihood
# i=1
for(i in 1:Nmc){
  Npnlike[i,] = dbinom(dat.2018$color.clean,1,x.ip_LOD[i,dat.2018$week-30],log=TRUE)  ## calculate log likelihoods
  Npnlike[i,is.na(Npnlike[i,])] = 0       ## missing data as weight 1; log(1)=0
  tmp = tapply(Npnlike[i,],dat.2018$week-30,sum)
  tmp2 = rep(0, length(unique(dat.2018$week)))
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
dat.2018_order = dat.2018[order(dat.2018$week),]
data_2018_loess_10 <- loess(as.numeric(as.character(color.clean)) ~ as.numeric(as.character(week)), data=dat.2018_order, span=0.10) # 10% smoothing span
data_2018_loess_10 <- predict(data_2018_loess_10) 
data_2018_loess_30 <- loess(as.numeric(as.character(color.clean)) ~ as.numeric(as.character(week)), data=dat.2018_order, span=0.30) # 10% smoothing span
data_2018_loess_30 <- predict(data_2018_loess_30) 

time <- 31:47
#Initial plotting of data and uncertainty partitioning BEFORE iteration is added
plot( dat.2018$week, dat.2018$color.clean ,pch="+",cex=0.5, xlab = "Day of Year", ylab = "Fall Color", main = "2018 Non-resampling Particle Filter")
ecoforecastR::ciEnvelope(time,ip_ci_LOD[1,],ip_ci_LOD[3,],col=col.alpha("purple",0.5))
ecoforecastR::ciEnvelope(time,i_ci_LOD[1,],i_ci_LOD[3,],col=col.alpha("green",0.5))
ecoforecastR::ciEnvelope(time,det_ci_LOD[1,],det_ci_LOD[3,],col=col.alpha("blue",0.5))
ecoforecastR::ciEnvelope(time,Npnpf[1,],Npnpf[3,],col=col.alpha("red",0.5))
lines(data_2018_loess_10, x=dat.2018_order$week, col="green", lwd = 2)
lines(data_2018_loess_30, x=dat.2018_order$week, col="blue", lwd = 2)


s = seq(31, 47, by = 2) # date that we're running the forecast on 
# set weights after the date to 0

iter_pf = list()
i=1
for(i in seq_along(s)) {
  # what row index corresponds to date for forecasting up to
  iter_pf[[i]] = matrix(NA,3,nobs)
  wbar = apply(Npnclike,2,mean)             ## mean weight at each time point
  tmp_clike = Npnclike
  tmp_clike[,(s[i]:47)-30] <- tmp_clike[,s[i]-30]
  for(j in 1:nobs){
    iter_pf[[i]][,j] = wtd.quantile(x.ip_LOD[,j],tmp_clike[,j]/wbar[j],c(0.025,0.5,0.975))  ## calculate weighted median and CI
  }
}

time <- 31:47
png(width= 750, filename= file.path(path.fig, paste0('Weekly Iterative fall color prediction', '.png')))
plot(dat.2018$week, dat.2018$color.clean ,pch="+",cex=0.5, xlab = "Week", ylab = "Fall Color", main = "Iterative 2018 Weekly Oak collection forecast")
ecoforecastR::ciEnvelope(time,ip_ci_LOD[1,],ip_ci_LOD[3,],col=col.alpha("purple",0.5))
for(i in seq_along(s)) {
  ecoforecastR::ciEnvelope(time[(s[i]:47)-30],iter_pf[[i]][1,(s[i]:47)-30],iter_pf[[i]][3,(s[i]:47)-30],col=col.alpha(i,0.5))
}
# ecoforecastR::ciEnvelope(time,Npnpf[1,],Npnpf[3,],col=col.alpha("red",0.5))
lines(data_2018_loess_10, x=dat.2018_order$week, col="green", lwd = 2)
lines(data_2018_loess_30, x=dat.2018_order$week, col="blue", lwd = 2)
dev.off()

