
#laod packages
library(jagsUI) #have to install jags directly on your machine
library(MCMCglmm)
require(plyr) #for ddply
require(tidyr) #for separate

#set working directory
# path.local <- "C:/Users/tosam/Documents/0_OSU/Dissertation/Data/SmallMammalCams/"
path.local <- "C:/Users/tosam/Documents/0_OSU/Dissertation/Manuscripts/SmallMammalCameras/SmammalCam_pub/"
#path.local <- "/raid1/home/fw/tosam/SmallMammalCams/"
setwd(path.local)

#write model file
cat("
      model
      {
      for (g in 1:n.sess) #for each grid
      {
      #### model for marked individuals
      
      for (i in 1:mmax) #individual i
      {
      zm[g,i] ~ dbern(psim[g]) 
      S[g,i,1] ~ dunif(xlims[1], xlims[2]) #activity center x coord
      S[g,i,2] ~ dunif(ylims[1], ylims[2]) #activity center y coord
      
      for(j in 1:J) #location j
      {
      D2[g,i,j] <- (S[g,i,1]-X[j,1])^2 + (S[g,i,2]-X[j,2])^2 #pythagorian theorem, distance to activity center
      lam[g,i,j] <- lam0[g]*exp(-D2[g,i,j]/(2*sigma[g]^2)) #lambda accounting for distance from activity center
      y[g,i,j] ~ dpois(lam.effm[g,i,j]*K[g,j])  #model accumulated counts across K for marked individuals and add for whether a camera is functioning or not
      lam.effm[g,i,j] <- lam[g,i,j]*zm[g,i] #lambda accounting for existence of individual
      }
      }
      #### priors for each grid
      log(sigma[g]) <- log.sigma[g] # take log sigma, note: now you cannot provide initial values for sigma
      log.sigma[g] ~ dnorm(mu.log.sigma, tau.log.sigma) # random effect on sigma for each grid
      
      log(lam0[g]) <- log.lam0[g] # take log lam0, note: now you cannot provide initial values for lam0
      log.lam0[g] ~ dnorm(mu.log.lam0, tau.log.lam0) #random effect on lam0 for each grid
      
      psim[g]~dunif(0,1)
      }
    
      for (g in 1:n.sess)
      {
      N[g] <- sum(zm[g,1:mmax])
      D[g] <- N[g]/S.Area
      }
    
      #### priors across grids
      #use log normal so all values are positive
      mu.log.sigma ~ dnorm(log(40), 1/(log(100)^2)) #mean sigma across all grids on log scale
      allgrids.mu.sigma <- exp(mu.log.sigma) #put on real scale
      sd.log.sigma ~ dunif(0,100) # sd of sigma on log (???) scale
      tau.log.sigma <- pow(sd.log.sigma,-2) # precision of sigma on log (???) scale
     
      mu.log.lam0 ~ dnorm(log(0.1), 1/(log(0.1)^2)) # mean lam0 across all grids on log scale (here, mean=0.1, sd=0.1 on real scale)
      allgrids.mu.lam <- exp(mu.log.lam0) #put on real scale
      sd.log.lam0 ~ dunif(0,5) # sd of lam0 on log (???) scale
      tau.log.lam0 <- pow(sd.log.lam0,-2) # precision of lam0 on log (???) scale

      } #end model description
      
      ", fill = TRUE, file="Models/SCR_multigrid_randomeffect.txt")

p <- function(sp.cr, n.chains, n.adapt, n.iter, n.burnin, n.cores)
{
  load(file=paste("RData/CR-SCR/SCR_randomeff_jagsargs_", sp.cr, ".RData", sep=""))
  
  mod.SCR_B <- jags(data=data.SCR_B, inits=inits, parameters.to.save = params, model.file="Models/SCR_multigrid_randomeffect.txt", n.chains=n.chains, 
                    n.adapt=n.adapt, n.iter=n.iter, n.burnin=n.burnin, n.thin=1, parallel = T, n.cores = n.cores)
  D.mode<-posterior.mode(as.mcmc(mod.SCR_B[[1]]$D))
  N.mode<-posterior.mode(as.mcmc(mod.SCR_B[[1]]$N))
  
  SCR_B <- list(sp.cr=sp.cr, mod=mod.SCR_B, Dhat=D.mode, Nhat=N.mode, SArea=data.SCR_B$S.Area)
  return(SCR_B)
}

pema <- p(sp.cr="PEMA", n.chains=3, n.adapt = 500, n.iter=2500, n.burnin=500, n.cores=3)
save(pema, file=paste(path.local, "RData/SCR_random_", pema$sp.cr, ".RData", sep=""))

glsa <- p(sp.cr="GLSA", n.chains=3, n.adapt = 500, n.iter=2500, n.burnin=500, n.cores=3)
save(glsa, file=paste(path.local, "RData/SCR_random_", glsa$sp.cr, ".RData", sep=""))

tato <- p(sp.cr="TATO", n.chains=3, n.adapt = 500, n.iter=2500, n.burnin=500, n.cores=3)
save(tato, file=paste(path.local, "RData/SCR_random_", tato$sp.cr, ".RData", sep=""))

###########

output <- pema
output <- glsa
output <- tato

out <- capture.output(output)
cat(out, file=paste("Output_SCR/", output$sp.cr, "_SCR-B_random_", Sys.Date(), ".txt", sep=""), sep="\n", append=F)

###########
#plot diagnostic plots

require(MCMCvis)

MCMCtrace(output$mod, params=c("lam0","sigma"), open_pdf=F)

#plot mcmc chains
jagsout.mcmc <- as.mcmc.list(output$mod$samples)
p <- xyplot(jagsout.mcmc[,c("allgrids.mu.lam", "allgrids.mu.sigma","N[1]","N[2]","N[3]","N[4]","N[5]","N[6]","N[7]","N[8]","N[9]")]) #plots traceplots in one figure
dp <- densityplot(jagsout.mcmc[,c("allgrids.mu.lam", "allgrids.mu.sigma","N[1]","N[2]","N[3]","N[4]","N[5]","N[6]","N[7]","N[8]","N[9]")]) #density plot of posterior

tiff(filename=paste("Figures_SCR/SCR-B_random_trace_", sp.cr, "_random_", Sys.Date(), ".tiff", sep=""), height=6, width=6, res=300, units="in", compression="lzw")
print(p)
dev.off()

tiff(filename=paste("Figures_SCR/SCR-B_random_density_", sp.cr, "_random_", Sys.Date(), ".tiff", sep=""), height=6, width=6, res=300, units="in", compression="lzw")
print(dp)
dev.off()

################
require(MCMCvis)

############
#load jags output
load(paste(path.local, "RData/SCR_random_TATO.RData", sep=""))

pema$SArea
# 5.022983

tato$SArea
# 5.022983

glsa$SArea
# 55.87544

s.all <- NULL
s2.all <- NULL
estimates <- function(model)
{
  s <- MCMCsummary(model$mod, params=c("lam0","sigma","psim","N","D"), round=3)
  s$param <- row.names(s)
  s$model <- "SCR_random"
  s$sp <- model$sp.cr
  
  s.all <<- rbind(s.all, s)
  
  MCMCtrace(model$mod, params=c("allgrids.mu.lam","allgrids.mu.sigma","psim","N"), type="density", filename=paste("Figures_SCR/SCR-B_density_random_", model$sp.cr,".pdf", sep=""))
  MCMCtrace(model$mod, params=c("allgrids.mu.lam","allgrids.mu.sigma","psim","N"), type="trace", filename=paste("Figures_SCR/SCR-B_trace_random_", model$sp.cr,".pdf", sep=""))
  
  s2 <- data.frame(dhat=model$Dhat, nhat=model$Nhat, model="SCR_random", sp=model$sp.cr, grid=1:length(model$Dhat))
  s2.all <<- rbind(s2.all, s2)
}

estimates(pema)
estimates(glsa)
estimates(tato)

# write.table(s.all, file="Output_SCR/SCR-B_random.txt", sep=",")
# write.table(s2.all, file="Output_SCR/SCR-B_random_mode.txt", sep=",")

s.all <- read.table("Output_SCR/SCR-B_random.txt", sep=",")
s2.all <- read.table("Output_SCR/SCR-B_random_mode.txt", sep=",")

scrrandom <- cbind(s.all[grep(s.all$param, pattern="D"),], s2.all[,c("dhat","nhat","grid")])
scrrandom$range <- paste(round(scrrandom$X2.5., digits=2), "-", round(scrrandom$X97.5., digits=2), sep=" ")
write.csv(scrrandom, file="../SmallMammalCams/Data_Final_ModelEstimates/SCR_random.csv")

bigpar <- NULL
pars <- function(model)
{
  s <- MCMCsummary(model$mod, params=c("allgrids.mu.lam","allgrids.mu.sigma"), round=3)
  s$param <- row.names(s)
  s$model <- "SCR_random"
  s$sp <- model$sp.cr
  
  bigpar <<- rbind(bigpar, s)
}
pars(pema)
pars(glsa)
pars(tato)

bigpar
#     sp      model mean.lam sd.lam 2.5%.lam 50%.lam 97.5%.lam Rhat.lam n.eff.lam       param.lam mean.sig sd.sig 2.5%.sig 50%.sig 97.5%.sig Rhat.sig n.eff.sig         param.sig
# 1 GLSA SCR_random    0.062  0.007    0.050   0.061     0.076     1.03        79 allgrids.mu.lam   39.460  1.826   36.084  39.372    42.985     1.04        55 allgrids.mu.sigma
# 2 PEMA SCR_random    0.063  0.013    0.042   0.062     0.090     1.00      3367 allgrids.mu.lam   12.834  1.409   10.244  12.769    15.822     1.00      6000 allgrids.mu.sigma
# 3 TATO SCR_random    0.168  0.022    0.127   0.166     0.213     1.00      2017 allgrids.mu.lam   33.522  1.438   30.755  33.500    36.493     1.00      6000 allgrids.mu.sigma

bigpar <- merge(bigpar[bigpar$param=="allgrids.mu.lam",], bigpar[bigpar$param=="allgrids.mu.sigma",], by=c("sp","model"), suffixes=c(".lam",".sig"))
# write.table(bigpar, file="C:/Users/tosam/Documents/0_OSU/Dissertation/Data/SmallMammalCams_SCR/SCR_bigparameters.txt", sep=",", row.names=F)
