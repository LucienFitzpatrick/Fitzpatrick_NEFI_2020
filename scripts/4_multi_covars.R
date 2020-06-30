library(coda)
library(rjags)
library(suncalc)

#########################
### load data, clean up #
#########################

path.doc <- ("../data_processed/fall/")

dat.npn <- read.csv(file.path(path.doc, "Fall_Phenology_data.csv"))


#creating 2018 frame for hindcasting
dat.2018 <- dat.npn[dat.npn$year == 2018, ]


#Setting the start of possible fall color as starting August 1st
dat.npn <- dat.npn[dat.npn$day_of_year > 213, ]

#Make 0,1 response variable numeric
dat.npn$color.full <- as.numeric(as.character(dat.npn$color.full))

# Calc CDD
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

# Calc length of day
lat <- 41.81405
long <- -88.05032
TZ_name <- "America/Chicago"
dayLengths <- numeric()
# days = unique(dat.npn$observation_date)
days = seq(as.Date("2019-08-02"), as.Date("2019-12-31"), by="days")
for(d in days){
  dte <- as.Date(d,origin="2019-08-01")
  suntimes <- getSunlightTimes(date=dte,
                               lat=lat,lon=long,keep=c("nauticalDawn","nauticalDusk"),
                               tz = TZ_name)
  dayLengths <- c(dayLengths,as.numeric(suntimes$nauticalDusk-suntimes$nauticalDawn))
}

time <- 214:365

covar_binom = "
model{

#### Data Model
for(i in 1:n){
y[i] ~ dbern(x[time[i]])
}

#### Process Model
for(t in 2:nt){
z[t]~dnorm(mu[t],tau_add)
mu[t] <- x[t-1]  + betaLOD*LOD[t] + betaCDD*CDD[t] # consider: betaX*x[t-1] + betaIntercept 
x[t] <- min(0.999,max(x[t-1],z[t]))
}

#### Priors
x[1] ~ dlnorm(x_ic,tau_ic)
# tau_obs ~ dgamma(a_obs,r_obs)
tau_add ~ dgamma(a_add,r_add)
betaCDD ~ dnorm(0, 1)
betaLOD ~ dnorm(0, 100)
}
"

day <- time-213
LOD.2019 <- as.data.frame(cbind(day, dayLengths))
CDD.2019 <- as.data.frame(cbind(day, CDD.2019))


data <- list(y = dat.npn$color.full, n = length(dat.npn$color.full), time = dat.npn$day_of_year-213, nt = 365-213, 
             a_add=0.1, r_add=0.001, x_ic = -100, tau_ic = 1000)


data$LOD = LOD.2019$dayLengths[match(data$time,LOD.2019$day)]
data$CDD = CDD.2019$CDD.2019[match(data$time,CDD.2019$day)]

#CDD_mod_out <- readRDS("model_output/CDDModel_Output.RDS")
#summary(CDD_mod_out)

# nchain = 3
# init <- list()
# for(i in 1:nchain){
#   init[[i]] <- list(betaCDD=rnorm(1,-0.004633,0.002013),tau_add=rnorm(1, 909.8, 0.00201))
# } 


j.model.covar   <- jags.model (file = textConnection(covar_binom),
                              data = data,
                              #inits = init,
                              n.chains = 3)

jags.out.covar   <- coda.samples (model = j.model.covar,
                                 variable.names = c("x","tau_add", "betaCDD", "betaLOD"),
                                 n.iter = 60000)
gelman.diag(jags.out.covar)
GBR <- gelman.plot(jags.out.covar)
burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[,,2]>1.1,1,any)),1)+1]


#gelman.plot(jags.out.covar[,"betacovar"])
#plot(jags.out.covar[,"x[151]"])


jags.burn.covar <- window(jags.out.covar,start=10000)  ## remove burn-in
gelman.diag(jags.burn.covars)                             ## check diagnostics post burn-in

#jags.out2 <- coda.samples(j.model,variable.names = c("b","S"),10000) # if you need additional samples

# acfplot(jags.burn.covar)
# 
effectiveSize(jags.burn.covar)
# 
cumuplot(jags.burn.covar,probs=c(0.025,0.25,0.5,0.75,0.975))

dic.samples(j.model.covar, 15000)


data_2019 = dat.npn[,c("day_of_year","color.full")]
data_2019 = data_2019[order(data_2019$day_of_year),]
data_2019$color.full = as.numeric(as.character(data_2019$color.full))
data_2019$day_of_year = as.numeric(as.character(data_2019$day_of_year))
data_2019_loess_10 <- loess(as.numeric(as.character(color.full)) ~ as.numeric(as.character(day_of_year)), data=data_2019, span=0.10) # 10% smoothing span
data_2019_loess_10 <- predict(data_2019_loess_10) 
data_2019_loess_30 <- loess(as.numeric(as.character(color.full)) ~ as.numeric(as.character(day_of_year)), data=data_2019, span=0.30) # 10% smoothing span
data_2019_loess_30 <- predict(data_2019_loess_30) 


out.covar <- as.matrix(jags.burn.covar)
x.cols <- grep("^x",colnames(out.covar)) ## grab all columns that start with the letter x
ci.covar <- apply(out.covar[,x.cols],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

time <- 214:365
plot(time,ci.covar[2,],type='n',ylim=c(0,1),ylab="Fall color")
## adjust x-axis label to be monthly if zoomed
ecoforecastR::ciEnvelope(time,ci.covar[1,],ci.covar[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(dat.npn$day_of_year, dat.npn$color.full ,pch="+",cex=0.5)
lines(data_2019_loess_30, x=data_2019$day_of_year, col="blue", lwd = 2)
