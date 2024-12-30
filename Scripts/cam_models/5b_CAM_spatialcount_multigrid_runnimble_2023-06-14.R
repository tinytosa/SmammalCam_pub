
#Spatial count models for camera trap data with different consolidation windows t=0, 15, 60, 1440

##############################
#load packages
require(nimble)
require(MCMCvis)

##############################
#load data

#import lambda0 and sigma values from SCR-B models, SCR-random
pars <- read.table(file="C:/Users/tosam/Documents/0_OSU/Dissertation/Data/SmallMammalCams_SCR/SCR_bigparameters.txt", sep=",", header=T)
pars <- pars[,c("sp","model","mean.lam","sd.lam","mean.sig","sd.sig")]

###############
#R defines dgamma as: f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s)
#mean=a*s and variance=a*s^2.
#mu = a*s --> a = mu/s
#v = a*s^2 --> v = mu*s = sd^2 --> s = sd^2/mu
#a = shape, s = scale

s <- function(mu, sd) {sd^2/mu}
a <- function(mu, s) {mu/s}

##########
#For Nimble
noninformative_multigrid_code <- nimbleCode({
  #priors across grids
  lam0~dunif(0,5)
  sigma~dunif(0,200)
  
  for (g in 1:n.sess) #for each grid
  {
    #priors for each grid
    psi[g]~dunif(0,1)
    
    ####model for unmarked individuals
    for (i in 1:M) #for each individual in each grid g
    {
      z[g,i]~dbern(psi[g]) 
      Su[g,i,1]~dunif(xlims[1], xlims[2])
      Su[g,i,2]~dunif(ylims[1], ylims[2])
      
      for(j in 1:J) #for each location j for each individual i in each grid g
      {
        D2u[g,i,j] <- (Su[g,i,1]-X[j,1])^2 + (Su[g,i,2]-X[j,2])^2
        lamu[g,i,j] <- lam0 * exp(-D2u[g,i,j]/(2*sigma^2)) # * z[g,i] * K[g,j]
        
        #yu[i,j]~dpois(lam.eff[i,j])
        lam.eff[g,i,j] <- lamu[g,i,j] * z[g,i] * K[g,j] # *Eff[j,k]  #add in whether camera was functioning or not
      }
    }
    
    for (j in 1:J) #for each trap in each grid
    {
      bigLambda[g,j] <- sum(lam.eff[g, 1:M, j])
      n[g,j] ~ dpois(lambda=bigLambda[g,j])
    }
  }
  
  for(g in 1:n.sess)
  {
    N[g] <- sum(z[g, 1:M]) 
    D[g] <- N[g]/S.Area*10000
  }
})

informative_multigrid_code <- nimbleCode({
  #priors
  lam0~dgamma(l.a, 1/l.s)
  sigma~dgamma(s.a, 1/s.s)
  
  for (g in 1:n.sess) #for each grid
  {
    #priors for each grid
    psi[g]~dunif(0,1)
    ####model for unmarked individuals
    for (i in 1:M) #for each individual in each grid g
    {
      z[g,i]~dbern(psi[g]) 
      Su[g,i,1]~dunif(xlims[1], xlims[2])
      Su[g,i,2]~dunif(ylims[1], ylims[2])
      
      for(j in 1:J) #for each location j for each individual i in each grid g
      {
        D2u[g,i,j]<-(Su[g,i,1]-X[j,1])^2 + (Su[g,i,2]-X[j,2])^2
        lamu[g,i,j]<-lam0*exp(-D2u[g,i,j]/(2*sigma^2)) # * z[g,i] * K[g,j]
        
        #yu[i,j]~dpois(lam.eff[i,j])
        lam.eff[g,i,j]<-lamu[g,i,j]*z[g,i]*K[g,j] # * Eff[j,k]  #add in whether camera was functioning or not
      }
    }
    
    for (j in 1:J)
    {
      bigLambda[g,j] <- sum(lam.eff[g,1:M,j])
      n[g,j] ~ dpois(lambda=bigLambda[g,j])
    }
  }
  
  for(g in 1:n.sess)
  {
    N[g] <- sum(z[g, 1:M]) 
    D[g] <- N[g]/S.Area*10000
  }
})

##############################
#load Rdata
##############################

sp.cr <- "PEMA"
n.chains = 3
n.adapt = 100
n.iter = 750
n.burnin = 100
# t = 1440

p <- function(t, sp.cr, n.chains, n.adapt, n.iter, n.burnin, n.cores, info)
{
  set.seed(123)
  # load("RData/CAM-SC_multigrid_nimbleargs_PEMA_test.RData")
  load(paste("RData/CAM-SC_multigrid_nimbleargs_", sp.cr, "_t", t, ".RData", sep=""))
  
  if(info == "noinfo")
  {
    constants <<- list(n.sess=nrow(data.UNM$n), M=M, J=J, S.Area=S.Area, K=K) #for non-informative priors
    locs <- cbind(cam$ssDF[[1]]$X, cam$ssDF[[1]]$Y)
    mod.UNM <- nimbleMCMC(code=noninformative_multigrid_code,
                          constants=constants, data=data.UNM, monitors=monitors, inits=inits,
                          niter=n.iter, nburnin=n.burnin, nchains=n.chains, 
                          check=T, progressBar=T, samplesAsCodaMCMC=T, summary=F)
  }
  
  if(info == "info")
  {
    print(pars[pars$sp == sp.cr,]$mean.lam)
    print(pars[pars$sp == sp.cr,]$sd.lam)
    
    l.s <- s(pars[pars$sp == sp.cr,]$mean.lam, pars[pars$sp == sp.cr,]$sd.lam)
    l.a <- a(pars[pars$sp == sp.cr,]$mean.lam, pars[pars$sp == sp.cr,]$sd.lam)
    
    print(l.s)
    print(l.a)
    
    #
    print(pars[pars$sp == sp.cr,]$mean.sig)
    print(pars[pars$sp == sp.cr,]$sd.sig)
    
    s.s <- s(pars[pars$sp == sp.cr,]$mean.sig, pars[pars$sp == sp.cr,]$sd.sig)
    s.a <- a(pars[pars$sp == sp.cr,]$mean.sig, pars[pars$sp == sp.cr,]$sd.sig)
    
    print(s.s)
    print(s.a)
    
    constants <<- list(n.sess=nrow(data.UNM$n), M=M, J=J, S.Area=S.Area, K=K, l.s=l.s, l.a=l.a, s.s=s.s, s.a=s.a) #for informative priors
    locs <- cbind(cam$ssDF[[1]]$X, cam$ssDF[[1]]$Y)
    mod.UNM <- nimbleMCMC(code=informative_multigrid_code, 
                          constants=constants, data=data.UNM, monitors=monitors, inits=inits,
                          niter=n.iter, nburnin=n.burnin, nchains=n.chains, 
                          check=T, progressBar=T, samplesAsCodaMCMC=T, summary=F)
  }
  out.UNM <- list(mod=mod.UNM, SArea=S.Area, sp.cr=sp.cr, M=M)
  return(out.UNM)
}

#########
#no info
for(t in c(0,15,60,1440))
{
  print(t)
  print(Sys.time())
  # pema <- p(t=t, sp.cr="PEMA", n.chains=3, n.adapt = 100, n.iter=1000, n.burnin=100, n.cores=3, info="noinfo")
  pema <- p(t=t, sp.cr="PEMA", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500, n.cores=3, info="noinfo")
  print(Sys.time())
  save(pema, file=paste("RData/CAM-SC_multigrid_noinfo_PEMA_t", t,".RData", sep=""))
  
}

for(t in c(0,15,60,1440))
{
  print(t)
  print(Sys.time())
  tato <- p(t=t, sp.cr="TATO", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500, n.cores=3, info="noinfo")
  print(Sys.time())
  save(tato, file=paste("RData/CAM-SC_multigrid_noinfo_TATO_t",t,".RData", sep=""))
}

for(t in c(0,15,60,1440))
{
  print(t)
  print(Sys.time())
  glsa <- p(t=t, sp.cr="GLSA", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500, n.cores=3, info="noinfo")
  print(Sys.time())
  save(glsa, file=paste("RData/CAM-SC_multigrid_noinfo_GLSA_t",t,".RData", sep=""))
}


#########
#info
for(t in c(0,15,60,1440))
{
  print(t)
  print(Sys.time())
  pema <- p(t=t, sp.cr="PEMA", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500, n.cores=3, info="info")
  print(Sys.time())
  save(pema, file=paste("RData/CAM-SC_multigrid_info_PEMA_t", t,".RData", sep=""))
  
}

for(t in c(0,15,60,1440))
{
  print(t)
  print(Sys.time())
  tato <- p(t=t, sp.cr="TATO", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500, n.cores=3, info="info")
  print(Sys.time())
  save(tato, file=paste("RData/CAM-SC_multigrid_info_TATO_t",t,".RData", sep=""))
}

for(t in c(0,15,60,1440))
{
  print(t)
  print(Sys.time())
  glsa <- p(t=t, sp.cr="GLSA", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500, n.cores=3, info="info")
  print(Sys.time())
  save(glsa, file=paste("RData/CAM-SC_multigrid_info_GLSA_t",t,".RData", sep=""))
}
