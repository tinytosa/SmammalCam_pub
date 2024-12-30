
#need to do this the first time running nimble
#make sure to reinstall RTools for R version, then reinstall nimble package

# path <- Sys.getenv('PATH')
# newPath <- paste("C:\\Rtools\\bin;C:\\Rtools\\mingw_64\\bin;",
#                  path, sep = "")
# Sys.setenv(PATH = newPath)

#####################
#laod packages
library(jagsUI) #have to install jags directly on your machine
library(MCMCglmm)
require(MCMCvis) #for MCMCsummary
require(plyr) #for ddply
require(tidyr) #for separate
require(nimble)

#####################
#load data
####
#import lambda0 and sigma values from SCR-B models, SCR-random
pars <- read.table(file="C:/Users/tosam/Documents/0_OSU/Dissertation/Data/SmallMammalCams_SCR/SCR_bigparameters.txt", sep=",", header=T)
pars <- pars[,c("sp","model","mean.lam","sd.lam","mean.sig","sd.sig")]

pars <- pars[rep(seq_len(nrow(pars)), each = 9), ]
pars$g <- rep(1:9, 3)

###############
#R defines dgamma as: f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s)
#mean=a*s and variance=a*s^2.
#mu = a*s --> a = mu/s
#v = a*s^2 --> v = mu*s = sd^2 --> s = sd^2/mu
#a = shape, s = scale

s <- function(mu, sd) {sd^2/mu}
a <- function(mu, s) {mu/s}

###################
#define model file
###################

##########
#For Nimble
noninformative_code <- nimbleCode({
  #priors
  lam0~dunif(0,5)
  sigma~dunif(0,200)
  psi~dunif(0,1)
  
  ####model for unmarked individuals
  for (i in 1:M)
  {
    z[i]~dbern(psi) 
    Su[i,1]~dunif(xlims[1], xlims[2])
    Su[i,2]~dunif(ylims[1], ylims[2])
    
    for(j in 1:J)
    {
      D2u[i,j]<-(Su[i,1]-X[j,1])^2 + (Su[i,2]-X[j,2])^2
      lamu[i,j]<-lam0*exp(-D2u[i,j]/(2*sigma^2))
      
      #yu[i,j]~dpois(lam.eff[i,j])
      lam.eff[i,j]<-lamu[i,j]*z[i]*K[j] # *Eff[j,k]  #add in whether camera was functioning or not
    }
  }
  
  for (j in 1:J)
  {
    bigLambda[j] <- sum(lam.eff[1:M,j])
    n[j] ~ dpois(bigLambda[j])
  }
  
  N <- sum(z[1:M]) 
  D <- N/S.Area*10000
})

informative_code <- nimbleCode({
  #priors
  lam0~dgamma(l.a, 1/l.s)
  sigma~dgamma(s.a, 1/s.s)
  psi~dunif(0,1)
  
  ####model for unmarked individuals
  for (i in 1:M)
  {
    z[i]~dbern(psi)
    Su[i,1]~dunif(xlims[1], xlims[2])
    Su[i,2]~dunif(ylims[1], ylims[2])
    
    for(j in 1:J)
    {
      D2u[i,j]<-(Su[i,1]-X[j,1])^2 + (Su[i,2]-X[j,2])^2
      lamu[i,j]<-lam0*exp(-D2u[i,j]/(2*sigma^2))
      
      #yu[i,j]~dpois(lam.eff[i,j])
      lam.eff[i,j]<-lamu[i,j]*z[i]*K[j] # *Eff[j,k]  #add in whether camera was functioning or not; z = indicator of real or fake, K = number of occasions
    }
  }
  
  for (j in 1:J)
  {
    bigLambda[j] <- sum(lam.eff[1:M,j])
    n[j] ~ dpois(bigLambda[j])
  }
  
  N <- sum(z[1:M])
  D <- N/S.Area*10000
})

# sp.cr <- "PEMA"
# g <- 1
# n.chains=3
# n.adapt = 100
# n.iter=1000
# n.burnin=500

p <- function(g, sp.cr, n.chains, n.adapt, n.iter, n.burnin, info)
{
  load(file=paste("RData/CR-SC_separate_nimbleargs_", sp.cr, "_SM-0", g,".RData", sep=""))

  if(info == "noinfo")
  {
    constants <<- list(M=M, J=J, S.Area=S.Area) #for non-informative priors
    print(constants)
    mod.UNM <- nimbleMCMC(code=noninformative_code, constants=constants, data=data.UNM, monitors=monitors, inits=inits, #non-informative priors
                          niter=n.iter, nburnin=n.burnin, nchains=n.chains, progressBar=T, check=F, samplesAsCodaMCMC=T, summary=F)
  }
  
  if(info == "info")
  {
    print(pars[pars$sp == sp.cr & pars$g == g,]$mean.lam)
    print(pars[pars$sp == sp.cr & pars$g == g,]$sd.lam)
    
    l.s <- s(pars[pars$sp == sp.cr & pars$g == g,]$mean.lam, pars[pars$sp == sp.cr & pars$g == g,]$sd.lam)
    l.a <- a(pars[pars$sp == sp.cr & pars$g == g,]$mean.lam, pars[pars$sp == sp.cr & pars$g == g,]$sd.lam)
    
    print(l.s)
    print(l.a)
    
    #
    print(pars[pars$sp == sp.cr & pars$g == g,]$mean.sig)
    print(pars[pars$sp == sp.cr & pars$g == g,]$sd.sig)
    
    s.s <- s(pars[pars$sp == sp.cr & pars$g == g,]$mean.sig, pars[pars$sp == sp.cr & pars$g == g,]$sd.sig)
    s.a <- a(pars[pars$sp == sp.cr & pars$g == g,]$mean.sig, pars[pars$sp == sp.cr & pars$g == g,]$sd.sig)
    
    print(s.s)
    print(s.a)
    constants <<- list(M=M, J=J, S.Area=S.Area, l.s=l.s, l.a=l.a, s.s=s.s, s.a=s.a) #for informative priors
    
    print(constants)
    mod.UNM <- nimbleMCMC(code=informative_code, constants=constants, data=data.UNM, monitors=monitors, inits=inits, #informative priors
                          niter=n.iter, nburnin=n.burnin, nchains=n.chains, progressBar=T, check=F, samplesAsCodaMCMC=T, summary=F)
    
  }
  out.UNM <- list(mod=mod.UNM, SArea=S.Area, sp.cr=sp.cr, g=g, info=info, M=M)
  return(out.UNM)
}

#################
#separate noinfo
#ran on 3/7/2023 w/ M=1000
Sys.time()
for(g in c(1:5,7:9))
{
  print(g)
  assign(paste("pema_",g, sep=""), p(g=g, sp.cr="PEMA", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500, info="noinfo"))
  # save.image(paste("RData/CR-SC_info_separate_PEMA_SM-0", g, ".RData", sep="")),
  save(list=ls(pattern="pema"), file=paste("RData/CR-SC_separate_noinfo_PEMA_SM-0", g, ".RData", sep=""))
  print(Sys.time())
}
save(list=ls(pattern="pema"), file=paste("RData/CR-SC_separate_noinfo_PEMA.RData", sep=""))

#################
#separate info
#################
#ran on 3/7/2023 w/ M=1000
#need to rerun with M=2000
Sys.time()
for(g in c(1:5,7:9))
{
  print(g)
  assign(paste("pema_",g, sep=""), p(g=g, sp.cr="PEMA", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500, info="info"))
  # save.image(paste("RData/CR-SC_info_separate_PEMA_SM-0", g, ".RData", sep=""))
  save(list=ls(pattern="pema"), file=paste("RData/CR-SC_separate_info_PEMA_SM-0", g, ".RData", sep=""))
  print(Sys.time())
}

save(list=ls(pattern="pema"), file=paste("RData/CR-SC_separate_info_PEMA.RData", sep=""))

# load("RData/CR-SC_info_separate_PEMA_SM-09.RData")
# save(pema_1,pema_2,pema_3,pema_4,pema_5,pema_7,pema_8,pema_9,
#      file=paste("RData/CR-SC_info_separate_PEMA.RData", sep=""))

# load("RData/CR-SC_info_separate_PEMA.RData")

# ran on 3/8/2023 w/ M=1000
Sys.time()
for(g in c(1:9))
{
  print(g)
  assign(paste("glsa_",g, sep=""), p(g, sp.cr="GLSA", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500, info="info"))
  # save.image(paste("RData/CR-SC_info_separate_GLSA_SM-0", g, ".RData", sep=""))
  save(list=ls(pattern="glsa"), file=paste("RData/CR-SC_separate_info_GLSA_SM-0", g, ".RData", sep=""))
  # save(glsa_9, file=paste("RData/CR-SC_separate_GLSA_SM-0", g, ".RData", sep=""))
  print(Sys.time())
}

save(list=ls(pattern="glsa"), file=paste("RData/CR-SC_separate_info_GLSA.RData", sep=""))

# ran on 3/8/2023 w/ M=1000
#ran grids 1:8 with M=1000
Sys.time()
for(g in c(1:9))
{
  print(g)
  assign(paste("tato_",g, sep=""), p(g, sp.cr="TATO", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500, info="info"))
  # save.image(paste("RData/CR-SC_info_separate_TATO_SM-0", g, ".RData", sep=""))
  save(list=ls(pattern="tato"), file=paste("RData/CR-SC_separate_info_TATO_SM-0", g, ".RData", sep=""))
  # save(tato_9, file=paste("RData/CR-SC_separate_TATO_SM-0", g, ".RData", sep=""))
  print(Sys.time())
}

save(list=ls(pattern="tato"), file=paste("RData/CR-SC_separate_info_TATO.RData", sep=""))

# #rerun tato grid 8 for both info and no info (change K from 7 to 6)
# #run on 3/7/2022 @13:52
# load("RData/CR-SC_separate_TATO.RData")
# print(Sys.time())
# g <- 8
# assign(paste("tato_",g, sep=""), p(g, sp.cr="TATO", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500))
# save.image(paste("RData/CR-SC_separate_TATO_SM-0", g, ".RData", sep=""))
# print(Sys.time())
# save(tato_1,tato_2,tato_3,tato_4,tato_5,tato_6,tato_7,tato_8,tato_9,
#      file="RData/CR-SC_separate_TATO.RData")
# 
# #run on 3/7/2023 @12:33
# load("RData/CR-SC_info_separate_TATO.RData")
# print(Sys.time())
# g <- 8
# assign(paste("tato_",g, sep=""), p(g, sp.cr="TATO", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500))
# save.image(paste("RData/CR-SC_info_separate_TATO_SM-0", g, ".RData", sep=""))
# print(Sys.time())
# save(tato_1,tato_2,tato_3,tato_4,tato_5,tato_6,tato_7,tato_8,tato_9,
#      file="RData/CR-SC_info_separate_TATO.RData")
# 
# 
# # load("RData/CR-SC_info_separate_GLSA_SM-09.RData")
# save(glsa_1,glsa_2,glsa_3,glsa_4,glsa_5,glsa_6,glsa_7,glsa_8,glsa_9,
#      file="RData/CR-SC_info_separate_GLSA.RData")
# save(tato_1,tato_2,tato_3,tato_4,tato_5,tato_6,tato_7,tato_8,tato_9,
#      file="RData/CR-SC_info_separate_TATO.RData")