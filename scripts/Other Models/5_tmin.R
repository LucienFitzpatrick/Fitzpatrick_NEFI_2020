#########################
##### load libraries  ###
#########################
library(coda)
library(rjags)

#########################
### load data, clean up #
#########################

spath.doc <- ("../data_processed/fall/")

dat.npn <- read.csv(file.path(path.doc, "Fall_Phenology_data.csv"))

#creating 2018 frame for hindcasting
dat.2018 <- dat.npn[dat.npn$year == 2018, ]

#Make 0,1 response variable numeric
dat.npn$color.full <- as.numeric(as.character(dat.npn$color.full))

time <- 214:365

##########################
### add tmin to process ###
##########################

tmin_binom = "
model{

#### Data Model
for(i in 1:n){
y[i] ~ dbern(x[time[i]])
}

#### Process Model
for(t in 2:nt){
z[t]~dnorm(mu[t],tau_add)
mu[t] <- x[t-1]  + betaTmin*tmin[t] # consider: betaX*x[t-1] + betaIntercept 
x[t] <- min(0.999,max(x[t-1],z[t]))
}

#### Priors
x[1] ~ dlnorm(x_ic,tau_ic)
# tau_obs ~ dgamma(a_obs,r_obs)
tau_add ~ dgamma(a_add,r_add)
betaTmin ~ dnorm(0, 1000)
}
"

day <- time-213
min.temp <- as.data.frame(cbind(day, min.temp.dat))

data <- list(y = dat.npn$color.full, n = length(dat.npn$color.full), time = dat.npn$day_of_year-212, nt = 365-213, 
             a_add=0.1, r_add=0.001, x_ic = -100, tau_ic = 1000)

data$tmin = min.temp$min.temp.dat[match(data$time,min.temp$day)]


j.model.tmin   <- jags.model (file = textConnection(tmin_binom),
                         data = data,
                         n.chains = 3)

jags.out.tmin   <- coda.samples (model = j.model.tmin,
                            variable.names = c("x","tau_add", "betaTmin"),
                            n.iter = 10000)


GBR <- gelman.plot(jags.out.tmin)
burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[,,2]>1.1,1,any)),1)+1]

gelman.plot(jags.out.tmin[,"betaTmin"])
plot(jags.out.tmin[,"x[151]"])

                    
jags.burn.tmin <- window(jags.out.tmin,start=burnin)  ## remove burn-in
gelman.diag(jags.burn.tmin)                             ## check diagnostics post burn-in

#jags.out2 <- coda.samples(j.model,variable.names = c("b","S"),10000) # if you need additional samples

# acfplot(jags.burn.tmin)
# 
effectiveSize(jags.burn.tmin)
# 
cumuplot(jags.burn.tmin,probs=c(0.025,0.25,0.5,0.75,0.975))

dic.samples(j.model.tmin, 15000)


data_2019 = dat.npn[,c("day_of_year","color.full")]
data_2019 = data_2019[order(data_2019$day_of_year),]
data_2019$color.full = as.numeric(as.character(data_2019$color.full))
data_2019$day_of_year = as.numeric(as.character(data_2019$day_of_year))
data_2019_loess_10 <- loess(as.numeric(as.character(color.full)) ~ as.numeric(as.character(day_of_year)), data=data_2019, span=0.10) # 10% smoothing span
data_2019_loess_10 <- predict(data_2019_loess_10) 
data_2019_loess_30 <- loess(as.numeric(as.character(color.full)) ~ as.numeric(as.character(day_of_year)), data=data_2019, span=0.30) # 10% smoothing span
data_2019_loess_30 <- predict(data_2019_loess_30) 


out.tmin <- as.matrix(jags.burn.tmin)
x.cols <- grep("^x",colnames(out.tmin)) ## grab all columns that start with the letter x
ci.tmin <- apply(out.tmin[,x.cols],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

time <- 214:365
plot(time,ci.tmin[2,],type='n',ylim=c(0,1),ylab="Fall color")
## adjust x-axis label to be monthly if zoomed
ecoforecastR::ciEnvelope(time,ci.tmin[1,],ci.tmin[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(dat.npn$day_of_year, dat.npn$color.full ,pch="+",cex=0.5)
lines(data_2019_loess_30, x=data_2019$day_of_year, col="blue", lwd = 2)

tmin.norm <- (min.temp.dat-min(min.temp.dat))/(max(min.temp.dat)-min(min.temp.dat))
#points(time,tmin.norm, ylab="Minimum temperature (C)")

saveRDS(out.tmin, "model_output/Tmin_Out_Matrix.RDS")
saveRDS(jags.burn.tmin, "model_output/Tmin_Out_JAGS.RDS")
