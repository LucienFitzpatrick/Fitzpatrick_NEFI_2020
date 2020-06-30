
path.doc <- ("../data_processed/fall/")

dat.npn <- read.csv(file.path(path.doc, "Fall_Phenology_data.csv"))

#creating 2018 frame for hindcasting
dat.2018 <- dat.npn[dat.npn$year == 2018, ]

dat.npn$color.full <- as.numeric(as.character(dat.npn$color.full))


Precip.2019 <- dat.met[dat.met$year==2019,]
Precip.2019<- Precip.2019[214:365,4]

time <- 214:365
day <- time-213
Precip.2019 <- as.data.frame(cbind(day, Precip.2019))


plot(Precip.2019$day, Precip.2019$Precip.2019, ylab="Precip")


RandomWalk_binom_Precip = "
model{

#### Data Model
for(i in 1:n){
y[i] ~ dbern(x[time[i]])
}

#### Process Model
for(t in 2:nt){
z[t]~dnorm(mu[t],tau_add)
mu[t] <- x[t-1]  + betaPrecip*Precip[t] # consider: betaX*x[t-1] + betaIntercept 
x[t] <- min(0.999,max(x[t-1],z[t]))
}

#### Priors
x[1] ~ dlnorm(x_ic,tau_ic)
# tau_obs ~ dgamma(a_obs,r_obs)
tau_add ~ dgamma(a_add,r_add)
betaPrecip ~ dnorm(0, .001)
}
"

data <- list(y = dat.npn$color.full, n = length(dat.npn$color.full), time = dat.npn$day_of_year-213, nt = 365-213, 
             a_add=1, r_add=0.00001, x_ic = -100, tau_ic = 1000)

data$Precip = Precip.2019$Precip.2019[match(data$time,Precip.2019$day)]


precip.model   <- jags.model (file = textConnection(RandomWalk_binom_Precip),
                         data = data,
                         n.chains = 3)

precip.out   <- coda.samples (model = precip.model,
                            variable.names = c("x", "tau_add", "betaPrecip"),
                            n.iter = 10000)

DIC.Precip <- dic.samples(precip.model, 5000)

gelman.diag(precip.out)

#TO DO: burn in!
burnin = 5000                                ## determine convergence
precip.burn <- window(precip.out,start=burnin)  ## remove burn-in
plot(precip.burn)                             ## check diagnostics post burn-in

p.out <- as.matrix(precip.burn)
x.cols <- grep("^x",colnames(p.out)) ## grab all columns that start with the letter x
precip.ci <- apply(p.out[,x.cols],2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

time <- 214:365
plot(time,precip.ci[2,],type='n',ylim=c(0,1),ylab="Fall color")
## adjust x-axis label to be monthly if zoomed
ecoforecastR::ciEnvelope(time, precip.ci[1,], precip.ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(dat.npn$day_of_year, dat.npn$color.full ,pch="+",cex=0.5)

saveRDS(precip.burn,"C:/Users/lucie/Documents/GitHub/NEFI/model_output/Precipmodel_Output.RDS")
saveRDS(DIC.Precip,"C:/Users/lucie/Documents/GitHub/NEFI/model_output/Precipmodel_DIC.RDS")


plot(p.out[,"tau_add"])




