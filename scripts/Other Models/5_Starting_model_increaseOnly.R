#########################
##### load libraries  ###
#########################
library(coda)
library(rjags)
library(ggplot2)

#########################
### load data, clean up #
#########################
path.doc <- ("../data_processed/fall/")
dir.create("../data_processed/model_output/", recursive = T, showWarnings = F)

path.doc <- ("../data_processed/fall/")

dat.npn <- read.csv(file.path(path.doc, "Fall_Phenology_data.csv"))

#creating 2018 frame for hindcasting
dat.2018 <- dat.npn[dat.npn$year == 2018, ]

#Setting the start of possible fall color as starting August 1st
dat.npn <- dat.npn[dat.npn$day_of_year > 213, ]

#Make 0,1 response variable numeric
dat.npn$color.full <- as.numeric(as.character(dat.npn$color.full))

#########################
#### Loess Smoothing ####
#########################

## 2019
data_2019 = dat.npn[,c("day_of_year","color.full")]
data_2019 = data_2019[order(data_2019$day_of_year),]
data_2019$color.full = as.numeric(as.character(data_2019$color.full))
data_2019$day_of_year = as.numeric(as.character(data_2019$day_of_year))

data_2019_loess_10 <- loess(as.numeric(as.character(color.full)) ~ as.numeric(as.character(day_of_year)), data=data_2019, span=0.10) # 10% smoothing span
data_2019_loess_10 <- predict(data_2019_loess_10) 
data_2019_loess_30 <- loess(as.numeric(as.character(color.full)) ~ as.numeric(as.character(day_of_year)), data=data_2019, span=0.30) # 10% smoothing span
data_2019_loess_30 <- predict(data_2019_loess_30) 
plot(data_2019$day_of_year, data_2019$color.full, type="l", main="Loess Smoothing 2019 Data", xlab="Day of Year", ylab="Fall Color")

#add these lines to figures below...
lines(data_2019_loess_10, x=data_2019$day_of_year, col="red", lwd = 2)
lines(data_2019_loess_30, x=data_2019$day_of_year, col="blue", lwd = 2)

# Calculate chilling degree days...

min.temp <- dat.met[dat.met$year==2019,]
min.temp <- min.temp[214:365,8]

CDD.2019 <- 0
for(i in 2:length(min.temp)){
  if((min.temp)[i] < 10){
    offset <- 10- (min.temp)[i]
  }else{
    offset <- 0
  }
  CDD.2019 <- c(CDD.2019,(CDD.2019[(i-1)]+offset))
}

time <- 214:365
plot(time,CDD.2019, ylab="Cold Degree Days")

##############################################
####### First try: Random walk        ########
##############################################

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

RW.data <- list(y = dat.npn$color.full, n = length(dat.npn$color.full), time = dat.npn$day_of_year-213, nt = 365-213, 
                a_obs=1, r_obs=0.0001, a_add=1, r_add=0.00001, x_ic = 0 , tau_ic = 1000)

RW.model   <- jags.model (file = textConnection(RandomWalk),
                          data = RW.data,
                          n.chains = 3)

RW.out   <- coda.samples (model = RW.model,
                          variable.names = c("x","tau_add"),
                          n.iter = 10000)

saveRDS(RW.out, "../data_processed/model_output/RandomWalk_Output.RDS")

RW.burnin = 3000                                ## determine convergence
RW.burn <- window(RW.out,start=RW.burnin)  ## remove burn-in

RW.DIC <- dic.samples(RW.model, 5000)
saveRDS(RW.DIC, "../data_processed/model_output/RandomWalk_DIC.RDS")

rw.out <- as.matrix(RW.out)
x.cols <- grep("^x",colnames(rw.out)) ## grab all columns that start with the letter x
rw.ci <- apply(rw.out[,x.cols],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

time.rng = c(1,length(RW.data$time))

jpeg("figures/Data.jpg")
plot(time,rw.ci[2,],type='n',ylim=c(0,1),ylab="Fall color")
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
}
points(dat.npn$day_of_year, dat.npn$color.full ,pch="+",cex=0.5)
lines(data_2019_loess_30, x=data_2019$day_of_year, col="blue", lwd = 2)
dev.off()


jpeg("figures/RandomWalk.jpg")
plot(time,rw.ci[2,],type='n',ylim=c(0,1),ylab="Fall color")
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
}
ecoforecastR::ciEnvelope(time,rw.ci[1,],rw.ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(dat.npn$day_of_year, dat.npn$color.full ,pch="+",cex=0.5)
lines(data_2019_loess_30, x=data_2019$day_of_year, col="blue", lwd = 2)
dev.off()


########################################
####### Random walk binomial    ########
########################################

RandomWalk_binom = "
model{

#### Data Model
for(i in 1:n){
y[i] ~ dbern(x[time[i]])
}

#### Process Model
for(t in 2:nt){
z[t]~dnorm(x[t-1],tau_add)
x[t] <- min(0.999,max(0.0001,z[t]))
}

#### Priors
x[1] ~ dnorm(x_ic,tau_ic)
# tau_obs ~ dgamma(a_obs,r_obs)
tau_add ~ dgamma(a_add,r_add)
}
"

RWB.data <- list(y = dat.npn$color.full, n = length(dat.npn$color.full), time = dat.npn$day_of_year-213, nt = 365-213, 
             a_add=1, r_add=0.00001, x_ic = 0 , tau_ic = 1000)

RWB.model   <- jags.model (file = textConnection(RandomWalk_binom),
                         data = RWB.data,
                         n.chains = 3)

#RWB.out   <- coda.samples (model = RWB.model,
#                            variable.names = c("x","tau_add"),
#                            n.iter = 10000)

#saveRDS(RWB.out, "../data_processed/model_output/RandomWalkBinom_Output.RDS")
RWB.out <- readRDS("../data_processed/model_output/RandomWalkBinom_Output.RDS")

RWB.burnin = 3000                                ## determine convergence
RWB.burn <- window(RWB.out,start=RWB.burnin)  ## remove burn-in

#RWB.DIC <- dic.samples(RWB.model, 5000)
#saveRDS(RWB.DIC, "../data_processed/model_output/RandomWalkBinom_DIC.RDS")
RWB.DIC <- readRDS("../data_processed/model_output/RandomWalkBinom_DIC.RDS")

rwb.out <- as.matrix(RWB.out)
x.cols <- grep("^x",colnames(rwb.out)) ## grab all columns that start with the letter x
rwb.ci <- apply(rwb.out[,x.cols],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

time.rng = c(1,length(RWB.data$time))

jpeg("figures/RandomWalk_Binom.jpg")
plot(time,rwb.ci[2,],type='n',ylim=c(0,1),ylab="Fall color")
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
}
ecoforecastR::ciEnvelope(time,rwb.ci[1,],rwb.ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(dat.npn$day_of_year, dat.npn$color.full ,pch="+",cex=0.5)
lines(data_2019_loess_30, x=data_2019$day_of_year, col="blue", lwd = 2)
dev.off()

###########################################
######### second try: phenology doesn't ###
######### move backward!      #############
###########################################

# Calling this "FWD" because phenology moves forward!

FWD = "
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

FWD.data <- list(y = dat.npn$color.full, n = length(dat.npn$color.full), time = dat.npn$day_of_year-213, nt = 365-213, 
            a_add=1, r_add=0.00001, x_ic = -100, tau_ic = 1000)


FWD.model   <- jags.model (file = textConnection(FWD),
                         data = FWD.data,
                         n.chains = 3)

FWD.out   <- coda.samples (model = FWD.model,
                            variable.names = c("x","tau_add"),
                            n.iter = 10000)

saveRDS(FWD.out, "../data_processed/model_output/PhenologyForward_Output.RDS")
FWD.out <- readRDS("../data_processed/model_output/PhenologyForward_Output.RDS")

FWD.burnin = 3000                                ## determine convergence
FWD.burn <- window(FWD.out,start=FWD.burnin)  ## remove burn-in

# any convergence problem?
out = list(params=NULL,predict=NULL) #Split output into parameters and state variables
mfit = as.matrix(FWD.out,chains=TRUE)
pred.cols = grep("x[",colnames(mfit),fixed=TRUE)
chain.col = which(colnames(mfit)=="CHAIN")
out$params = ecoforecastR::mat2mcmc.list(mfit[,-pred.cols])
GBR.vals <- gelman.diag(out$params)
GBR.vals
# no.

#FWD.DIC <- dic.samples(FWD.model, 5000)
#saveRDS(FWD.DIC, "../data_processed/model_output/PhenologyForward_DIC.RDS")
FWD.DIC <- readRDS("../data_processed/model_output/PhenologyForward_DIC.RDS")

fwd.out <- as.matrix(FWD.out)
x.cols <- grep("^x",colnames(fwd.out)) ## grab all columns that start with the letter x
fwd.ci <- apply(fwd.out[,x.cols],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

jpeg("figures/Phenology_Forward.jpg")
plot(time,fwd.ci[2,],type='n',ylim=c(0,1),ylab="Fall color")
## adjust x-axis label to be monthly if zoomed
ecoforecastR::ciEnvelope(time,fwd.ci[1,],fwd.ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(dat.npn$day_of_year, dat.npn$color.full ,pch="+",cex=0.5)
lines(data_2019_loess_30, x=data_2019$day_of_year, col="blue", lwd = 2)
dev.off()


##########################
### add CDD to process ###
##########################

CDD = "
model{

#### Data Model
for(i in 1:n){
y[i] ~ dbern(x[time[i]])
}

#### Process Model
for(t in 2:nt){
z[t]~dnorm(mu[t],tau_add)
mu[t] <- x[t-1]  + betaCDD*CDD[t] # consider: betaX*x[t-1] + betaIntercept 
x[t] <- min(0.999,max(x[t-1],z[t]))
}

#### Priors
x[1] ~ dlnorm(x_ic,tau_ic)
# tau_obs ~ dgamma(a_obs,r_obs)
tau_add ~ dgamma(a_add,r_add)
betaCDD ~ dnorm(0, 1)
}
"

day <- time-213
CDD.2019 <- as.data.frame(cbind(day, CDD.2019))

CDD.data <- list(y = dat.npn$color.full, n = length(dat.npn$color.full), time = dat.npn$day_of_year-213, nt = 365-213, 
             a_add=1, r_add=0.00001, x_ic = -100, tau_ic = 1000)

CDD.data$CDD = CDD.2019$CDD.2019[match(CDD.data$time,CDD.2019$day)]


CDD.model   <- jags.model (file = textConnection(CDD),
                         data = CDD.data,
                         n.chains = 3)

#CDD.out   <- coda.samples (model = CDD.model,
#                            variable.names = c("x","tau_add", "betaCDD"),
#                            n.iter = 10000)

#saveRDS(CDD.out, "../data_processed/model_output/CDDModel_Output.RDS")
CDD.out <- readRDS("../data_processed/model_output/CDDModel_Output.RDS")

CDD.burnin = 3000                                ## determine convergence
CDD.burn <- window(CDD.out,start=CDD.burnin)  ## remove burn-in

#CDD.DIC <- dic.samples(CDD.model, 5000)
#saveRDS(CDD.DIC, "../data_processed/model_output/CDDModel_DIC.RDS")
CDD.DIC <- readRDS("../data_processed/model_output/CDDModel_DIC.RDS")


cdd.out <- as.matrix(CDD.out)
x.cols <- grep("^x",colnames(cdd.out)) ## grab all columns that start with the letter x
cdd.ci <- apply(cdd.out[,x.cols],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

jpeg("figures/CDDModel.jpg")
plot(time,cdd.ci[2,],type='n',ylim=c(0,1),ylab="Fall color")
## adjust x-axis label to be monthly if zoomed
ecoforecastR::ciEnvelope(time,cdd.ci[1,],cdd.ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(dat.npn$day_of_year, dat.npn$color.full ,pch="+",cex=0.5)
lines(data_2019_loess_30, x=data_2019$day_of_year, col="blue", lwd = 2)
dev.off()


################################
### add short-wave radiation ###
################################

SWR = "
model{

#### Data Model
for(i in 1:n){
y[i] ~ dbern(x[time[i]])
}

#### Process Model
for(t in 2:nt){
z[t]~dnorm(mu[t],tau_add)
mu[t] <- x[t-1]  + betaCDD*CDD[t] # consider: betaX*x[t-1] + betaIntercept 
x[t] <- min(0.999,max(x[t-1],z[t]))
}

#### Priors
x[1] ~ dlnorm(x_ic,tau_ic)
# tau_obs ~ dgamma(a_obs,r_obs)
tau_add ~ dgamma(a_add,r_add)
betaCDD ~ dnorm(0, 1)
}
"

day <- time-213
CDD.2019 <- as.data.frame(cbind(day, CDD.2019))

CDD.data <- list(y = dat.npn$color.full, n = length(dat.npn$color.full), time = dat.npn$day_of_year-213, nt = 365-213, 
                 a_add=1, r_add=0.00001, x_ic = -100, tau_ic = 1000)

CDD.data$CDD = CDD.2019$CDD.2019[match(CDD.data$time,CDD.2019$day)]


CDD.model   <- jags.model (file = textConnection(CDD),
                           data = CDD.data,
                           n.chains = 3)

CDD.out   <- coda.samples (model = CDD.model,
                           variable.names = c("x","tau_add", "betaCDD"),
                           n.iter = 10000)

saveRDS(CDD.out, "../data_processed/model_output/CDDModel_Output.RDS")

CDD.burnin = 3000                                ## determine convergence
CDD.burn <- window(CDD.out,start=CDD.burnin)  ## remove burn-in

CDD.DIC <- dic.samples(CDD.model, 5000)
saveRDS(CDD.DIC, "../data_processed/model_output/CDDModel_DIC.RDS")

cdd.out <- as.matrix(CDD.out)
x.cols <- grep("^x",colnames(cdd.out)) ## grab all columns that start with the letter x
cdd.ci <- apply(cdd.out[,x.cols],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

jpeg("figures/CDDModel.jpg")
plot(time,cdd.ci[2,],type='n',ylim=c(0,1),ylab="Fall color")
## adjust x-axis label to be monthly if zoomed
ecoforecastR::ciEnvelope(time,cdd.ci[1,],cdd.ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(dat.npn$day_of_year, dat.npn$color.full ,pch="+",cex=0.5)
lines(data_2019_loess_30, x=data_2019$day_of_year, col="blue", lwd = 2)
dev.off()

######################################
##### adding plots to one another ####
######################################
jpeg("figures/RW.jpg")
plot(time,rw.ci[2,],type='n',ylim=c(0,1),ylab="Fall color")
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
}
ecoforecastR::ciEnvelope(time,rwb.ci[1,],rwb.ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.9))
points(dat.npn$day_of_year, dat.npn$color.full ,pch="+",cex=0.5)
lines(data_2019_loess_30, x=data_2019$day_of_year, col="blue", lwd = 2)
dev.off()

#add phenology forward
jpeg("figures/RWplusFWD.jpg")
plot(time,rw.ci[2,],type='n',ylim=c(0,1),ylab="Fall color")
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
}
ecoforecastR::ciEnvelope(time,rwb.ci[1,],rwb.ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.9))
ecoforecastR::ciEnvelope(time,fwd.ci[1,],fwd.ci[3,],col=ecoforecastR::col.alpha("Maroon",0.75))
points(dat.npn$day_of_year, dat.npn$color.full ,pch="+",cex=0.5)
lines(data_2019_loess_30, x=data_2019$day_of_year, col="blue", lwd = 2)
dev.off()

#add CDD
jpeg("figures/FWDplusCDD.jpg")
plot(time,rw.ci[2,],type='n',ylim=c(0,1),ylab="Fall color")
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
}
ecoforecastR::ciEnvelope(time,fwd.ci[1,],fwd.ci[3,],col=ecoforecastR::col.alpha("Maroon",0.75))
ecoforecastR::ciEnvelope(time,cdd.ci[1,],cdd.ci[3,],col=ecoforecastR::col.alpha("Gold",0.5))
points(dat.npn$day_of_year, dat.npn$color.full ,pch="+",cex=0.5)
lines(data_2019_loess_30, x=data_2019$day_of_year, col="blue", lwd = 2)
dev.off()
