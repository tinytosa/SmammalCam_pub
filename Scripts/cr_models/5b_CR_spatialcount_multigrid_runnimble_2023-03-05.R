
##############################
#load packages
require(nimble)
require(MCMCvis)

##############################
#load data

#import lambda0 and sigma values from SCR-B models, SCR-random
pars <- read.table(file="C:/Users/tosam/Documents/0_OSU/Dissertation/Data/SmallMammalCams_SCR/SCR_bigparameters.txt", sep=",", header=T)
pars <- pars[,c("sp","model","mean.lam","sd.lam","mean.sig","sd.sig")]

# pars <- pars[rep(seq_len(nrow(pars)), each = 9), ]
# pars$g <- rep(1:9, 3)

###############
#R defines dgamma as: f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s)
#mean=a*s and variance=a*s^2.
#mu = a*s --> a = mu/s
#v = a*s^2 --> v = mu*s = sd^2 --> s = sd^2/mu
#a = shape, s = scale

s <- function(mu, sd) {sd^2/mu}
a <- function(mu, s) {mu/s}

##########
#Nimble models
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
        D2u[g,i,j]<-(Su[g,i,1]-X[j,1])^2 + (Su[g,i,2]-X[j,2])^2
        lamu[g,i,j]<-lam0*exp(-D2u[g,i,j]/(2*sigma^2))
        
        #yu[i,j]~dpois(lam.eff[i,j])
        lam.eff[g,i,j]<-lamu[g,i,j]*z[g,i]*K[g,j] # *Eff[j,k]  #add in whether camera was functioning or not
      }
    }
    
    for (j in 1:J)
    {
      bigLambda[g,j] <- sum(lam.eff[g,1:M,j])
      n[g,j] ~ dpois(bigLambda[g,j])
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
        lamu[g,i,j]<-lam0*exp(-D2u[g,i,j]/(2*sigma^2))
        
        #yu[i,j]~dpois(lam.eff[i,j])
        lam.eff[g,i,j]<-lamu[g,i,j]*z[g,i]*K[g,j] # *Eff[j,k]  #add in whether camera was functioning or not
      }
    }
    
    for (j in 1:J)
    {
      bigLambda[g,j] <- sum(lam.eff[g,1:M,j])
      n[g,j] ~ dpois(bigLambda[g,j])
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

# sp.cr <- "PEMA"
# n.chains=3
# n.adapt = 1000
# n.iter = 7500
# n.burnin = 1000

p <- function(sp.cr, n.chains, n.adapt, n.iter, n.burnin, n.cores, info)
{
  set.seed(123)
  load(paste("RData/CR-SC_nimbleargs/CR-SC_multigrid_nimbleargs_", sp.cr,".RData", sep=""))
  
  if(info == "noinfo")
  {
    constants <<- list(n.sess=nrow(data.UNM$n), M=M, J=J, S.Area=S.Area) #for non-informative priors
    mod.UNM <- nimbleMCMC(code=noninformative_multigrid_code,
                          constants=constants, data=data.UNM, monitors=monitors, inits=inits,
                          niter=n.iter, nburnin=n.burnin, nchains=n.chains, 
                          check=F, progressBar=T, samplesAsCodaMCMC=T, summary=F)
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
    
    constants <<- list(n.sess=nrow(data.UNM$n), M=M, J=J, S.Area=S.Area, l.s=l.s, l.a=l.a, s.s=s.s, s.a=s.a) #for informative priors
    mod.UNM <- nimbleMCMC(code=informative_multigrid_code, 
                          constants=constants, data=data.UNM, monitors=monitors, inits=inits,
                          niter=n.iter, nburnin=n.burnin, nchains=n.chains, 
                          check=F, progressBar=T, samplesAsCodaMCMC=T, summary=F)
  }
  out.UNM <- list(mod=mod.UNM, SArea=S.Area, sp.cr=sp.cr, M=M)
  return(out.UNM)
}

#no info
#ran on 3/7/2023
#reran on 3/9/2023 with M=1000 on chinook
Sys.time()
pema <- p(sp.cr="PEMA", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500, n.cores=3, info="noinfo")
Sys.time()
save(pema, file="RData/CR-SC_multigrid_noinfo_PEMA.RData")

#ran on 3/7/2023
Sys.time()
tato <- p(sp.cr="TATO", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500, n.cores=3, info="noinfo")
Sys.time()
save(tato, file="RData/CR-SC_multigrid_noinfo_TATO.RData")

#ran on 3/7/2023
Sys.time()
glsa <- p(sp.cr="GLSA", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500, n.cores=3, info="noinfo")
Sys.time()
save(glsa, file="RData/CR-SC_multigrid_noinfo_GLSA.RData")

####################
#info

#need to run
Sys.time()
pema <- p(sp.cr="PEMA", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500, n.cores=3, info="info")
Sys.time()
save(pema, file="RData/CR-SC_multigrid_info_PEMA.RData")

#need to run
Sys.time()
tato <- p(sp.cr="TATO", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500, n.cores=3, info="info")
Sys.time()
save(tato, file="RData/CR-SC_multigrid_info_TATO.RData")

#need to run
Sys.time()
glsa <- p(sp.cr="GLSA", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500, n.cores=3, info="info")
Sys.time()
save(glsa, file="RData/CR-SC_multigrid_info_GLSA.RData")

#####################
results <- as.mcmc.list(mod.UNM)
modes <- posterior.mode(results)
rhat <- gelman.diag(results, multivariate = F)
# posterior.mode(results$chain1)
# posterior.mode(results$chain2)
# posterior.mode(results$chain3)

out.UNM<-list(mod=mod.UNM, results=results, modes=modes, rhat=rhat, SArea=S.Area)
return(out.UNM)
