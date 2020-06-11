
# path.hub <- "C:/Users/lucie/Documents/GitHub/NEFI/data/"

path.hub <- "D:/git_proj/Fitzpatrick_NEFI_2020/data"

dat.npn <- read.csv(file.path(path.hub, file = "Arb_Quercus_NPN_data_leaves_CLEAN_individual.csv"), na.strings = "-9999")

#Daymet for when using covariates
dat.met <- read.csv(file.path(path.hub, file = "Daymet_data_raw.csv"))

#creating 2018 frame for hindcasting
dat.2018 <- dat.npn[dat.npn$year == 2018, ]

#isolating just 2019 year for our model
dat.npn <- dat.npn[dat.npn$year == 2019, ]

#Setting the start of possible fall color as starting August 1st
dat.npn <- dat.npn[dat.npn$day_of_year > 213, ]


dat.npn$color.full <- as.numeric(as.character(dat.npn$color.full))


library(rjags)
library(coda)


RandomWalk = "
model{
  
  #### Data Model
  for(i in 1:n){
    y[i] ~ dnorm(x[time[i]], tau_obs)

  }
  
  #### Process Model
  for(t in 2:nt){
    x[t]~dnorm(x[t-1],tau_add)
  }
  
  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add)
}
"


RandomWalk_binom = "
model{
  
  #### Data Model
  for(i in 1:n){
    y[i] ~ dbern(x[time[i]])
  }
  
  #### Process Model
  for(t in 2:nt){
    z[t]~dnorm(x[t-1],tau_add)
    x[t] <- min(0.999,max(x[t-1],z[t]))
  }
  
  #### Priors
  x[1] ~ dlnorm(x_ic,tau_ic)
  # tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add)
}
"
data <- list(y = dat.npn$color.full, n = length(dat.npn$color.full), time = dat.npn$day_of_year-213, nt = 365-213, 
             a_add=1, r_add=1, x_ic = -100 , tau_ic = 1000)




j.model   <- jags.model (file = textConnection(RandomWalk_binom),
                         data = data,
                         n.chains = 3)

jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add"),
                            n.iter = 10000)


gelman.diag(jags.out)

#plot(jags.out)

#GBR <- gelman.plot(jags.out)

#burnin = 5000                                ## determine convergence
#jags.burn <- window(jags.out,start=burnin)  ## remove burn-in
#plot(jags.burn)                             ## check diagnostics post burn-in

out <- as.matrix(jags.out)
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

time <- 214:365

plot.run <- function(ci) {
  plot(time,ci[2,],type='n',ylim=c(0,1),ylab="Fall color")
  ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
  points(dat.npn$day_of_year, dat.npn$color.full ,pch="+",cex=0.5)
}

#### Forecast function 
##` @param IC    Initial Conditions
##` @param Q     Process error (default = 0 for deterministic runs)
##` @param n     Size of Monte Carlo ensemble

NT = data$nt
forecastx <- function(IC,Q=0,n=Nmc){
  x <- z <- matrix(NA,n,NT)  ## storage
  Xprev <- IC           ## initialize
  for(t in 1:NT){
    z[,t] = rnorm(n,Xprev,Q)
    x[,t] <- pmin(0.999,pmax(Xprev,z[,t]))
    Xprev <- x[,t]                                  ## update IC
  }
  return(x)
}


## initial conditions
Nmc = 1000
IC <- rlnorm(Nmc, data$x_ic,(1/sqrt(data$tau_ic)))

x.det <- forecastx(IC=mean(IC),
                   Q=mean((1/sqrt(out[,"tau_add"])))/10,  ## process error off
                   n=Nmc)

x.i <- forecastx(IC=IC,
                 Q=mean((1/sqrt(out[,"tau_add"])))/10,  ## process error off
                 n=Nmc)

x.ip <- forecastx(IC=IC,
                  Q=1/sqrt(out[,"tau_add"])/10,  ## process error off
                  n=Nmc)


## Plot run
det_ci <- apply(x.det,2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
plot.run(ci)
lines(time,x.det,col="purple",lwd=3)

i_ci <- apply(x.i,2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
plot.run(ci)
lines(time,x.i,col="purple",lwd=3)

ip_ci <- apply(x.ip,2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

plot(dat.npn$day_of_year, dat.npn$color.full ,pch="+",cex=0.5, xlab = "Day of Year", ylab = "Fall Color", main = "2019")
ecoforecastR::ciEnvelope(time,ip_ci[1,],ip_ci[3,],col=col.alpha("purple",0.5))
ecoforecastR::ciEnvelope(time,i_ci[1,],i_ci[3,],col=col.alpha("green",0.5))
ecoforecastR::ciEnvelope(time,det_ci[1,],det_ci[3,],col=col.alpha("blue",0.5))

# 2018 data, not 2019
# scoring? 
# anything to see the trend in the data itself? (smoothing curve (loess), binned means)


#### Forecasting 2018
dat.2018
data_2018 <- list(y = dat.2018$color.full, n = length(dat.2018$color.full), time = dat.2018$day_of_year-213, nt = 365-213, 
             a_add=1, r_add=1, x_ic = -100 , tau_ic = 1000)

Nmc = 1000
IC_2018 <- rlnorm(Nmc, data_2018$x_ic,(1/sqrt(data_2018$tau_ic)))

x.det_2018 <- forecastx(IC=mean(IC),
                   Q=mean((1/sqrt(out[,"tau_add"])))/10,  ## process error off
                   n=Nmc)

x.i_2018 <- forecastx(IC=IC,
                 Q=mean((1/sqrt(out[,"tau_add"])))/10,  ## process error off
                 n=Nmc)

x.ip_2018 <- forecastx(IC=IC,
                  Q=1/sqrt(out[,"tau_add"])/10,  ## process error off
                  n=Nmc)

det_ci_2018 <- apply(x.det_2018,2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

i_ci_2018 <- apply(x.i_2018,2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

ip_ci_2018 <- apply(x.ip_2018,2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

plot(dat.npn$day_of_year, dat.npn$color.full ,pch="+",cex=0.5, xlab = "Day of Year", ylab = "Fall Color", main = "2018")
ecoforecastR::ciEnvelope(time,ip_ci_2018[1,],ip_ci_2018[3,],col=col.alpha("purple",0.5))
ecoforecastR::ciEnvelope(time,i_ci_2018[1,],i_ci_2018[3,],col=col.alpha("green",0.5))
ecoforecastR::ciEnvelope(time,det_ci_2018[1,],det_ci_2018[3,],col=col.alpha("blue",0.5))


#### Loess Smoothing
library(ggplot2)

## 2019
data_2019 = dat.npn[,c("day_of_year","color.full")]
data_2019 = data_2019[order(data_2019$day_of_year),]
data_2019$color.full = as.numeric(as.character(data_2019$color.full))
data_2019$day_of_year = as.numeric(as.character(data_2019$day_of_year))

plot(data_2019$day_of_year, data_2019$color.full, type = "l")

data_2019_loess_10 <- loess(as.numeric(as.character(color.full)) ~ as.numeric(as.character(day_of_year)), data=data_2019, span=0.10) # 10% smoothing span
data_2019_loess_10 <- predict(data_2019_loess_10) 
data_2019_loess_30 <- loess(as.numeric(as.character(color.full)) ~ as.numeric(as.character(day_of_year)), data=data_2019, span=0.30) # 10% smoothing span
data_2019_loess_30 <- predict(data_2019_loess_30) 

plot(data_2019$day_of_year, data_2019$color.full, type="l", main="Loess Smoothing 2019 Data", xlab="Day of Year", ylab="Fall Color")
lines(data_2019_loess_10, x=data_2019$day_of_year, col="red", lwd = 2)
lines(data_2019_loess_30, x=data_2019$day_of_year, col="blue", lwd = 2)

plot(data_2019_loess_10, x=data_2019$day_of_year, col="red", lwd = 2, type = 'l')
lines(data_2019_loess_30, x=data_2019$day_of_year, col="blue", lwd = 2)

## 2018
data_2018= dat.2018[,c("day_of_year","color.full")]
data_2018 = data_2018[order(data_2018$day_of_year),]
data_2018$color.full = as.numeric(as.character(data_2018$color.full)) # 4 NAs
which(is.na(data_2018$color.full))
data_2018[c(1,2,34,36),]
data_2018 = data_2018[-c(1,2,34,36),]
data_2018$day_of_year = as.numeric(as.character(data_2018$day_of_year))

plot(data_2018$day_of_year, data_2018$color.full, type = "l")

data_2018_loess_10 <- loess(as.numeric(as.character(color.full)) ~ as.numeric(as.character(day_of_year)), data=data_2018, span=0.10) # 10% smoothing span
data_2018_loess_10 <- predict(data_2018_loess_10) 
data_2018_loess_30 <- loess(as.numeric(as.character(color.full)) ~ as.numeric(as.character(day_of_year)), data=data_2018, span=0.30) # 10% smoothing span
data_2018_loess_30 <- predict(data_2018_loess_30) 

plot(data_2018$day_of_year, data_2018$color.full, type="l", main="Loess Smoothing 2018 Data", xlab="Day of Year", ylab="Fall Color")
lines(data_2018_loess_10, x=data_2018$day_of_year, col="red", lwd = 2)
lines(data_2018_loess_30, x=data_2018$day_of_year, col="blue", lwd = 2)

plot(data_2018_loess_10, x=data_2018$day_of_year, col="red", lwd = 2, type = 'l')
lines(data_2018_loess_30, x=data_2018$day_of_year, col="blue", lwd = 2)