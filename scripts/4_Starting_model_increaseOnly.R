dat.npn <- read.csv("Arb_Quercus_NPN_data_leaves_CLEAN_individual.csv", na.strings = "-9999")

#Daymet for when using covariates
dat.met <- read.csv("Daymet_data_raw.csv")

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
  x[1] ~ dlnorm(x_ic,tau_ic)
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
            a_add=1, r_add=1, x_ic = -100, tau_ic = 1000)




j.model   <- jags.model (file = textConnection(RandomWalk_binom),
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

#TO DO: burn in.
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

#######################
## chilling degrees  ##
#######################
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

plot(time,CDD.2019)

##########################
### add CDD to process ###
##########################

RandomWalk_binom = "
model{

#### Data Model
for(i in 1:n){
y[i] ~ dbern(x[time[i]])
}

#### Process Model
for(t in 2:nt){
z[t]~dnorm(x[t-1],tau_add)
mu[t] <- x[t-1]  + betaX*x[t-1] + betaIntercept*Xf[t,1] + betaTmin*Xf[t,2]
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

#TO DO: match dates for data$y and CDD.2019, make sure dimensions are the same across lists
#TO DO pt 2: Look at Exercise 06 together make sure we know where the covariates go
data$CDD <- CDD.2019

ef.out <- ecoforecastR::fit_dlm(model=list(obs="y",fixed="~ 1 + X + CDD"),data)



j.model   <- jags.model (file = textConnection(RandomWalk_binom),
                         data = data,
                         n.chains = 3)

jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add"),
                            n.iter = 10000)