

##########################
#spatial count model from capture recapture with small mammal camera trap data
#estimates for each grid calculated separately
#run SC model chains in parallel
##########################

######
#load packages

# library(oSCR)
# require(coda)
require(nimble)
# require(igraph)
# require(plyr)
require(MCMCglmm)
# require(lattice)

#######
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

noninformative_code <- nimbleCode({
  #priors
  lam0~dunif(0,100)
  sigma~dunif(0,200)
  psi~dunif(0,1)
  
  ####model for unmarked individuals
  for (i in 1:M)
  {
    z[i]~dbern(psi) #is individual real? 
    Su[i,1]~dunif(xlims[1], xlims[2]) #x-coord activity center
    Su[i,2]~dunif(ylims[1], ylims[2]) #y-coord activity center
    
    for(j in 1:J)
    {
      D2u[i,j]<-(Su[i,1]-X[j,1])^2 + (Su[i,2]-X[j,2])^2 #squared distance
      lamu[i,j]<-lam0*exp(-D2u[i,j]/(2*sigma^2)) * z[i] #capture prob at trap j
      lam.eff[i,j]<-lamu[i,j]*K[j] # *Eff[j,k]  #add in whether camera was functioning or not, K is number of occasions each camera was operational
    }
  }
  
  for (j in 1:J)
  {
    bigLambda[j] <- sum(lam.eff[1:M,j])
    n[j] ~ dpois(bigLambda[j]) #model for data
  }
  
  N <- sum(z[1:M]) 
  D <- N/S.Area*10000 #density in ha
})

informative_code <- nimbleCode({
  #priors
  lam0~dunif(0,100)
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
      lamu[i,j]<-lam0*exp(-D2u[i,j]/(2*sigma^2)) * z[i]
      
      #yu[i,j]~dpois(lam.eff[i,j])
      lam.eff[i,j]<-lamu[i,j]*K[j] # *Eff[j,k]  #add in whether camera was functioning or not
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

########
#function to run in parallel

run_all_code <- function(seed, code, data, constants, monitors, inits, n.adapt, n.iter, n.burnin)
{
  library(nimble)
  
  params<- c("N","D","lambda0","sigma","psi")
  
  myModel <- nimbleModel(code = code,
                         data = data,
                         constants = constants,
                         inits = inits)
  CmyModel <- compileNimble(myModel)
  
  myMCMC <- configureMCMC(CmyModel, monitors=monitors)
  BmyMCMC <- buildMCMC(myMCMC)
  CmyMCMC <- compileNimble(BmyMCMC)
  
  # results <- runMCMC(CmyMCMC, niter = 7500, nburnin=1000, setSeed = seed, summary=F, samplesAsCodaMCMC = T, progressBar = T)
  results <- runMCMC(CmyMCMC, niter = n.iter, nburnin=n.burnin, setSeed = seed, summary=F, samplesAsCodaMCMC = T, progressBar = T)
  
  return(results)
}

p <- function(g, sp.cr, t, n.chains, n.adapt, n.iter, n.burnin, info)
{
  # g <- 3
  # sp.cr <- "PEMA"
  # info <- "info"
  # t <- 0
  
  load(file=paste("RData/CAM-SC_separate_nimbleargs_", sp.cr,"_SM-0", g,"_t",t,".RData", sep=""))
  print(rbind(data.UNM$n, K=K))
  
  if(info == "noinfo")
  {
    constants <<- list(M=M, K=K, J=J, S.Area=S.Area) #for non-informative priors
    # print(constants)
    
    #######
    #run in parallel
    this_cluster <- parallel::makeCluster(3)
    
    chain_output <- parallel::parLapply(cl = this_cluster, X = 1:3, 
                                        fun = run_all_code, 
                                        code=noninformative_code, #non-informative priors
                                        data = data.UNM,
                                        constants=constants,monitors=monitors,inits=inits,
                                        n.adapt=n.adapt, n.iter=n.iter, n.burnin=n.burnin)
    parallel::stopCluster(this_cluster)
  }
    
  if(info == "info")
  {
    l.s <- s(pars[pars$sp == sp.cr & pars$g == g,]$mean.lam, pars[pars$sp == sp.cr & pars$g == g,]$sd.lam)
    l.a <- a(pars[pars$sp == sp.cr & pars$g == g,]$mean.lam, pars[pars$sp == sp.cr & pars$g == g,]$sd.lam)

    s.s <- s(pars[pars$sp == sp.cr & pars$g == g,]$mean.sig, pars[pars$sp == sp.cr & pars$g == g,]$sd.sig)
    s.a <- a(pars[pars$sp == sp.cr & pars$g == g,]$mean.sig, pars[pars$sp == sp.cr & pars$g == g,]$sd.sig)

    constants <<- list(M=M, K=K, J=J, S.Area=S.Area, l.s=l.s, l.a=l.a, s.s=s.s, s.a=s.a) #for informative priors
    print(constants)
    
    #######
    #run in parallel
    this_cluster <- parallel::makeCluster(3)
    
    chain_output <- parallel::parLapply(cl = this_cluster, X = 1:3, 
                                        fun = run_all_code, 
                                        code=informative_code, #informative priors
                                        data = data.UNM,
                                        constants=constants,monitors=monitors,inits=inits,
                                        n.adapt=n.adapt, n.iter=n.iter, n.burnin=n.burnin)
    parallel::stopCluster(this_cluster)
    # mod.UNM <- nimbleMCMC(code=informative_code, constants=constants, data=data.UNM, monitors=monitors, inits=inits, #informative priors
    #                       niter=n.iter, nburnin=n.burnin, nchains=n.chains, progressBar=T, check=F, samplesAsCodaMCMC=T, summary=F)
  }
  out.UNM <- list(mod=chain_output, SArea=S.Area, sp.cr=sp.cr, g=g, info=info, M=M, t=t, J=J, K=K)
  return(out.UNM)
}


###############
#separate noinfo
###############

# g <- 1
# sp.cr <- "TATO"
# info <- "noinfo"
# t <- 0

#
#takes roughly 1 hour per grid for PEMA
Sys.time()
# t <- 0
# t <- 15
# t <- 60
t <- 1440
for(g in c(1:8))
{
  print(g)
  assign(paste("pema_",g, "_t",t, sep=""), p(g=g, t=t, sp.cr="PEMA", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500, info="noinfo"))
  save(list=ls(pattern="pema"), file=paste("RData/CAM-SC_separate_noinfo_PEMA_SM-0", g,"_t", t, ".RData", sep=""))
  print(Sys.time())
}

save(list=ls(pattern="pema"), file=paste("RData/CAM-SC_separate_noinfo_PEMA_t", t, ".RData", sep=""))

#
#takes roughly 2 hours per grid for TATO
Sys.time()
# t <- 0
for(t in c(0,15,60,1440))
{
  print(t)
  for(g in c(1:8))
  {
    print(g)
    assign(paste("tato_",g, "_t",t, sep=""), p(g=g, t=t, sp.cr="TATO", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500, info="noinfo"))
    save(list=ls(pattern="tato"), file=paste("RData/CAM-SC_separate_noinfo_TATO_SM-0", g,"_t", t, ".RData", sep=""))
    print(Sys.time())
  }
  save(list=ls(pattern="tato"), file=paste("RData/CAM-SC_separate_noinfo_TATO_t", t, ".RData", sep=""))
}


#
#takes roughly 2.5 hours per grid for GLSA
Sys.time()
# t <- 0
for(t in c(0,15,60,1440))
{
  print(t)
  for(g in c(1:8))
  {
    print(g)
    assign(paste("glsa_",g, "_t",t, sep=""), p(g=g, t=t, sp.cr="GLSA", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500, info="noinfo"))
    save(list=ls(pattern="glsa"), file=paste("RData/CAM-SC_separate_noinfo_GLSA_SM-0", g,"_t", t, ".RData", sep=""))
    print(Sys.time())
  }
  save(list=ls(pattern="glsa"), file=paste("RData/CAM-SC_separate_noinfo_GLSA_t", t, ".RData", sep=""))
}

##############################
#separate info
##############################

#
Sys.time()
# t <- 0
# t <- 15
# t <- 60
# t <- 1440
for(g in c(1:8))
{
  print(g)
  assign(paste("pema_",g, "_t",t, sep=""), p(g=g, t=t, sp.cr="PEMA", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500, info="info"))
  save(list=ls(pattern="pema"), file=paste("RData/CAM-SC_separate_info_PEMA_SM-0", g,"_t", t, ".RData", sep=""))
  print(Sys.time())
}
save(list=ls(pattern="pema"), file=paste("RData/CAM-SC_separate_info_PEMA_t", t, ".RData", sep=""))

#
Sys.time()
# t <- 0
for(t in c(0,15,60,1440))
{
  print(t)
  for(g in c(1:8))
  {
    print(g)
    assign(paste("tato_",g, "_t",t, sep=""), p(g=g, t=t, sp.cr="TATO", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500, info="info"))
    save(list=ls(pattern="tato"), file=paste("RData/CAM-SC_separate_info_TATO_SM-0", g,"_t", t, ".RData", sep=""))
    print(Sys.time())
  }
  save(list=ls(pattern="tato"), file=paste("RData/CAM-SC_separate_info_TATO_t", t, ".RData", sep=""))
}

#
Sys.time()
# t <- 0
for(t in c(0,15,60,1440))
{
  print(t)
  for(g in c(1:8))
  {
    print(g)
    assign(paste("glsa_",g, "_t",t, sep=""), p(g=g, t=t, sp.cr="GLSA", n.chains=3, n.adapt = 2500, n.iter=10000, n.burnin=2500, info="info"))
    save(list=ls(pattern="glsa"), file=paste("RData/CAM-SC_separate_info_GLSA_SM-0", g,"_t", t, ".RData", sep=""))
    print(Sys.time())
  }
  save(list=ls(pattern="glsa"), file=paste("RData/CAM-SC_separate_info_GLSA_t", t, ".RData", sep=""))
}

# run_all_code <- function(seed, data, constants, monitors, inits)
# {
#   library(nimble)
#   #library(nimble, lib.loc = "/raid1/home/fw/ruprechj/R/powerpc64le-redhat-linux-gnu-library/3.6")
#   
#   params<- c("N","D","lambda0","sigma","psi")
#   
#   myModel <- nimbleModel(code = noninformative_code,
#                          data = data,
#                          constants = constants,
#                          inits = inits)
#   CmyModel <- compileNimble(myModel)
#   
#   myMCMC <- configureMCMC(CmyModel, monitors=monitors)
#   BmyMCMC <- buildMCMC(myMCMC)
#   CmyMCMC <- compileNimble(BmyMCMC)
#   
#   results <- runMCMC(CmyMCMC, niter = 7500, nburnin=1000, setSeed = seed, summary=F, samplesAsCodaMCMC = T, progressBar = T)
#   
#   return(results)
#   
# }
# 
# ########
# #set up data to run
# #specify these
# g <- 1
# for(g in c(1:8))
# {
#   print(g)
#   #sp.cam <- "DEER MOUSE"
#   #sp.cr <- "PEMA"
#   sp.cam <- "TOWNSENDS CHIPMUNK"
#   sp.cr <- "TATO"
#   #sp.cam <- "FLYING SQUIRREL"
#   #sp.cr <- "GLSA"
#   t <- 1440
#   M <- 2500
#   #traps <- traps.cam #trap locs for cameras, 8x8 tomahawk grid + 5x5 sherman grid
#   #traps <- traps.cam.small
#   traps <- traps.large
#   
#   #
#   print(g)
#   print(sp.cam)
#   print(t)
#   
#   data.cam <- cam.data[which(cam.data$Keywords==sp.cam & cam.data$grid==g),] #subset data for species and grid
#   data.cam <- data.cam[data.cam$Loc %in% traps$RC,] #subset so that only data from cameras in "traps"
#   
#   cam <- format.fun(data=data.cam, traps=traps, sp=sp.cam, grid=g, t=t)
#   
#   n <- cam$sf$caphist[[1]][,,1]
#   S.Area <- nrow(cam$ssDF[[1]])*cam$res^2 #this in in meters
#   K <- K.all[K.all$grid == g,]$occ
#   
#   print(n)
#   print(S.Area)
#   print(K)
#   
#   #For Nimble
#   set.seed(123)
#   data.UNM <- list(n=n, X=cam$sf$traps[[1]],
#                    xlims=c(min(cam$ssDF[[1]]$X), max(cam$ssDF[[1]]$X)),
#                    ylims=c(min(cam$ssDF[[1]]$Y), max(cam$ssDF[[1]]$Y)))
#   constants <- list(M=M, J=ncol(cam$sf$caphist[[1]]), S.Area=S.Area, K=K) #for non-informative priors
#   #constants <- list(M=M, J=ncol(data$sf$caphist[[1]]), S.Area=S.Area, K=K, s.s=s.s, s.a=s.a) #for informative priors
#   monitors <- c("lam0", "sigma", "psi", "N", "D") #same as params
#   inits <- list(lam0=runif(1), sigma=runif(1), psi=runif(1), z=rbinom(M,1,0.5),
#                 Su=cbind(runif(M, min(cam$ssDF[[1]]$X), max(cam$ssDF[[1]]$X)), #starting location of all M
#                          runif(M, min(cam$ssDF[[1]]$Y), max(cam$ssDF[[1]]$Y))))
#   
#   #######
#   #run in parallel
#   
#   this_cluster <- parallel::makeCluster(3)
#   
#   chain_output <- parallel::parLapply(cl = this_cluster, X = 1:3, 
#                                       fun = run_all_code, 
#                                       data = data.UNM,
#                                       constants=constants,
#                                       monitors=monitors,
#                                       inits=inits)
#   parallel::stopCluster(this_cluster)
#   
#   results <- as.mcmc.list(chain_output)
#   modes <- posterior.mode(results)
#   rhat <- gelman.diag(results, multivariate = F) #seed must be different for each run
#   
#   out.UNM<-list(mod=results, results=results, modes=modes, rhat=rhat, SArea=S.Area)
#   
#   out <- capture.output(summary(out.UNM$mod), out.UNM$modes, out.UNM$rhat, out.UNM$SArea)
#   
#   save(out.UNM, file=paste("Output_",run, "/Results_", sp.cr, "_SM-0", g, "_", Sys.Date(), "_t", t,"_m", M, ".RData", sep=""))
#   
#   cat(out, file=paste("Output_", run, "/UNM-B-cam_separate_", sp.cr, "_SM-0", g, "_", Sys.Date(), "_t", t,"_m", M, ".txt", sep=""), sep="\n", append=F)
#   
#   tp <- xyplot(results) #plots traceplots in one figure
#   dp <- densityplot(results) #density plot of posterior
#   
#   tiff(filename=paste("Figures_", run, "/outnimble_trace_", sp.cr, "_SM-0", g, "_", Sys.Date(), "_t", t,"_m", M,".tiff", sep=""), height=6, width=6, res=300, units="in", compression="lzw")
#   print(tp)
#   dev.off()
#   
#   tiff(filename=paste("Figures_", run, "/outnimble_densityplot_", sp.cr, "_SM-0", g, "_", Sys.Date(),  "_t", t,"_m",M, ".tiff", sep=""), height=6, width=6, res=300, units="in", compression="lzw")
#   print(dp)
#   dev.off()
# }

# results <- as.mcmc.list(chain_output)
# modes <- posterior.mode(results)
# rhat <- gelman.diag(results, multivariate = F) #seed must be different for each run
# 
# out.UNM<-list(mod=results, results=results, modes=modes, rhat=rhat, SArea=S.Area)
# 
# out <- capture.output(summary(out.UNM$mod), out.UNM$modes, out.UNM$rhat, out.UNM$SArea)

# MCMCtrace(model$mod, params=c("lam0","sigma","psi","N","D"), type="density", ind=T, filename=paste("Figures_CR-SC/density_", priors,"_", type, "_", model$sp.cr,"_", model$g, ".pdf", sep=""), open_pdf=T)
# MCMCtrace(model$mod, params=c("lam0","sigma","psi","N","D"), type="trace", ind=T, filename=paste("Figures_CR-SC/trace_", priors, "_", type,"_", model$sp.cr,"_", model$g,".pdf", sep=""), open_pdf=F)