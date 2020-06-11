#########################
##### load libraries  ###
#########################
library(coda)
library(rjags)

#########################
### load data, clean up #
#########################

#Morton arboretum color change data from 2018-2019, NPN protocols
dat.npn <- read.csv("Arb_Quercus_NPN_data_leaves_CLEAN_individual.csv", na.strings = "-9999")

#Daymet for when using covariates
dat.met <- read.csv("Daymet_data_raw.csv")

#creating 2018 frame for hindcasting
dat.2018 <- dat.npn[dat.npn$year == 2018, ]

#isolating just 2019 year for our model
dat.npn <- dat.npn[dat.npn$year == 2019, ]

#Setting the start of possible fall color as starting August 1st
dat.npn <- dat.npn[dat.npn$day_of_year > 213, ]

#Make 0,1 response variable numeric
dat.npn$color.full <- as.numeric(as.character(dat.npn$color.full))


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
####### First try: Random walk Binom! ########
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
# tau_obs ~ dgamma(a_obs,r_obs)
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
x[t] <- min(0.999,max(0.0001,z[t]))
}

#### Priors
x[1] ~ dnorm(x_ic,tau_ic)
# tau_obs ~ dgamma(a_obs,r_obs)
tau_add ~ dgamma(a_add,r_add)
}
"
data <- list(y = dat.npn$color.full, n = length(dat.npn$color.full), time = dat.npn$day_of_year-213, nt = 365-213, 
             a_add=1, r_add=1, x_ic = 0 , tau_ic = 1000)

j.model   <- jags.model (file = textConnection(RandomWalk_binom),
                         data = data,
                         n.chains = 3)

jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add"),
                            n.iter = 10000)


#gelman.diag(jags.out)
#plot(jags.out)
#GBR <- gelman.plot(jags.out)

#burnin = 5000                                ## determine convergence
#jags.burn <- window(jags.out,start=burnin)  ## remove burn-in
#plot(jags.burn)                             ## check diagnostics post burn-in

out <- as.matrix(jags.out)
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

time.rng = c(1,length(data$time))
plot(time,ci[2,],type='n',ylim=c(0,1),ylab="Fall color")
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
}
ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(dat.npn$day_of_year, dat.npn$color.full ,pch="+",cex=0.5)

###########################################
######### second try: phenology doesn't ###
######### move backward!      #############

RandomWalk_binom_prog = "
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
            a_add=1, r_add=1, x_ic = -100, tau_ic = 1000)




j.model   <- jags.model (file = textConnection(RandomWalk_binom_prog),
                         data = data,
                         n.chains = 3)

jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add"),
                            n.iter = 10000)


gelman.diag(jags.out)

# is this a convergence problem?
out = list(params=NULL,predict=NULL) #Split output into parameters and state variables
mfit = as.matrix(jags.out,chains=TRUE)
pred.cols = grep("x[",colnames(mfit),fixed=TRUE)
chain.col = which(colnames(mfit)=="CHAIN")
out$params = ecoforecastR::mat2mcmc.list(mfit[,-pred.cols])
GBR.vals <- gelman.diag(out$params)
GBR.vals
# no.

#plot(jags.out)
#GBR <- gelman.plot(jags.out)

#TO DO: burn in!
#burnin = 5000                                ## determine convergence
#jags.burn <- window(jags.out,start=burnin)  ## remove burn-in
#plot(jags.burn)                             ## check diagnostics post burn-in

out <- as.matrix(jags.out)
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

time <- 214:365
plot(time,ci[2,],type='n',ylim=c(0,1),ylab="Fall color")
## adjust x-axis label to be monthly if zoomed
ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(dat.npn$day_of_year, dat.npn$color.full ,pch="+",cex=0.5)


##########################
### add CDD to process ###
##########################

RandomWalk_binom_CDD = "
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

data <- list(y = dat.npn$color.full, n = length(dat.npn$color.full), time = dat.npn$day_of_year-213, nt = 365-213, 
             a_add=1, r_add=0.00001, x_ic = -100, tau_ic = 1000)

data$CDD = CDD.2019$CDD.2019[match(data$time,CDD.2019$day)]


j.model   <- jags.model (file = textConnection(RandomWalk_binom_CDD),
                         data = data,
                         n.chains = 3)

jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add", "betaCDD"),
                            n.iter = 10000)

dic.samples(j.model, 5000)


#plot(jags.out)
#GBR <- gelman.plot(jags.out)

#TO DO: burn in!
#burnin = 5000                                ## determine convergence
#jags.burn <- window(jags.out,start=burnin)  ## remove burn-in
#plot(jags.burn)                             ## check diagnostics post burn-in

out <- as.matrix(jags.out)
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

time <- 214:365
plot(time,ci[2,],type='n',ylim=c(0,1),ylab="Fall color")
## adjust x-axis label to be monthly if zoomed
ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(dat.npn$day_of_year, dat.npn$color.full ,pch="+",cex=0.5)

