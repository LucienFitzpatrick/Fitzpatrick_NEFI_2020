#### Length of day model

library(suncalc)

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
plot(time, dayLengths)

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
x[t] <- min(0.999,max(x[t-1],z[t]))
}

#### Priors
x[1] ~ dlnorm(x_ic,tau_ic)
# tau_obs ~ dgamma(a_obs,r_obs)
tau_add ~ dgamma(a_add,r_add)
betaLOD ~ dnorm(0, 1000)
}
"
day <- time-213
LOD.2019 <- as.data.frame(cbind(day, dayLengths))

data <- list(y = dat.npn$color.full, n = length(dat.npn$color.full), time = dat.npn$day_of_year-213, nt = 365-213, 
             a_add=0.1, r_add=0.001, x_ic = -100, tau_ic = 1000)


data$LOD = LOD.2019$dayLengths[match(data$time,LOD.2019$day)]


j.model_LOD   <- jags.model (file = textConnection(binom_LOD),
                         data = data,
                         n.chains = 3)

jags.out_LOD   <- coda.samples (model = j.model_LOD ,
                            variable.names = c("x","tau_add", "betaLOD"),
                            n.iter = 10000)

# GBR_LOD <- gelman.plot(jags.out_LOD)
gelman.diag(jags.out_LOD) # tau_add high
gelman.plot(jags.out_LOD[,"betaLOD"])

# rank order should be consistent, plot one that's the slowest to converge  

burnin_LOD = 4000                                ## determine convergence
jags.burn_LOD <- window(jags.out_LOD,start=burnin)  ## remove burn-in
#plot(jags.burn) 
summary(jags.burn_LOD)

out_LOD <- as.matrix(jags.burn_LOD)
saveRDS(out_LOD, "model_output/LengthOfDay_Output.RDS")

x.cols_LOD <- grep("^x",colnames(out_LOD)) ## grab all columns that start with the letter x
ci_LOD <- apply(out_LOD[,x.cols_LOD],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

time <- 214:365
plot(time,ci_LOD[2,],type='n',ylim=c(0,1),ylab="Fall color")
## adjust x-axis label to be monthly if zoomed
ecoforecastR::ciEnvelope(time,ci_LOD[1,],ci_LOD[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(dat.npn$day_of_year, dat.npn$color.full ,pch="+",cex=0.5)
lines(data_2019_loess_10, x=data_2019$day_of_year, col="red", lwd = 2)
lines(data_2019_loess_30, x=data_2019$day_of_year, col="blue", lwd = 2)

dic_LOD = dic.samples(j.model_LOD, 5000)
saveRDS(dic_LOD, "model_output/LengthOfDay_DIC.RDS")



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




Nmc = 1000
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

det_ci_LOD <- apply(x.det_LOD,2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
i_ci_LOD <- apply(x.i_LOD,2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
ip_ci_LOD <- apply(x.ip_LOD,2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale


plot(dat.npn$day_of_year, dat.npn$color.full ,pch="+",cex=0.5, xlab = "Day of Year", ylab = "Fall Color", main = "2019")
ecoforecastR::ciEnvelope(time,ip_ci_LOD[1,],ip_ci_LOD[3,],col=col.alpha("purple",0.5))
ecoforecastR::ciEnvelope(time,i_ci_LOD[1,],i_ci_LOD[3,],col=col.alpha("green",0.5))
ecoforecastR::ciEnvelope(time,det_ci_LOD[1,],det_ci_LOD[3,],col=col.alpha("blue",0.5))



# 2018 forecast
dat.2018
dat.2018 <- dat.2018[dat.2018$day_of_year > 213, ]
dat.2018$day_of_year <- as.numeric(as.character(dat.2018$day_of_year))
dat.2018$color.full <- as.numeric(as.character(dat.2018$color.full))

data_2018 <- list(y = dat.2018$color.full, n = length(dat.2018$color.full), time = dat.2018$day_of_year-213, nt = 365-213, 
             a_add=0.1, r_add=0.001, x_ic = -100, tau_ic = 1000)

Nmc = 10000
IC <- rlnorm(Nmc, data_2018$x_ic,(1/sqrt(data_2018$tau_ic)))

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

det_ci_LOD <- apply(x.det_LOD,2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
i_ci_LOD <- apply(x.i_LOD,2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale
ip_ci_LOD <- apply(x.ip_LOD,2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale


plot(dat.2018$day_of_year, dat.2018$color.full ,pch="+",cex=0.5, xlab = "Day of Year", ylab = "Fall Color", main = "2018 Forecast")
ecoforecastR::ciEnvelope(time,ip_ci_LOD[1,],ip_ci_LOD[3,],col=col.alpha("purple",0.5))
ecoforecastR::ciEnvelope(time,i_ci_LOD[1,],i_ci_LOD[3,],col=col.alpha("green",0.5))
ecoforecastR::ciEnvelope(time,det_ci_LOD[1,],det_ci_LOD[3,],col=col.alpha("blue",0.5))
lines(data_2018_loess_10, x=dat.2018_order$day_of_year, col="red", lwd = 2)
lines(data_2018_loess_30, x=dat.2018_order$day_of_year, col="blue", lwd = 2)




#####################################
# Particle Filter

## calculate the cumulative likelihoods
## to be used as PF weights
Npnlike = matrix(NA,nrow = Nmc, ncol = nrow(dat.2018))
Npnclike = matrix(0,nrow = Nmc, ncol = 365-213)
# i=1
for(i in 1:Nmc){
  Npnlike[i,] = dbinom(dat.2018$color.full,1,x.ip_LOD[i,dat.2018$day_of_year-213],log=TRUE)  ## calculate log likelihoods
  Npnlike[i,is.na(Npnlike[i,])] = 0       ## missing data as weight 1; log(1)=0
  # Npnlike[i,] = exp(cumsum(Npnlike[i,]))  ## convert to cumulative likelihood
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

