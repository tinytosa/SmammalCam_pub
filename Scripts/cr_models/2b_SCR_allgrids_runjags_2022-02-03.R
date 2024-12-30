
#format data/create arguments for JAGS using 2a_SCR_allgrids_jagsargs_2022-02-03.R


######
#install.packages("devtools")
#library(devtools)
#install_github("jaroyle/oSCR")

######
#library(oSCR)
library(jagsUI) #have to install jags directly on your machine
library(MCMCglmm)
require(plyr) #for ddply
require(tidyr) #for separate

path.local <- "C:/Users/tosam/Documents/0_OSU/Dissertation/Manuscripts/SmallMammalCameras/SmammalCam_pub/"
#path.local <- "/raid1/home/fw/tosam/SmallMammalCams/"
setwd(path.local)

sink("Models/SCR_multigrid.txt")
cat("
      model
      {
      for (g in 1:n.sess) #for each grid
      {
      #### priors for each grid
      psim[g]~dunif(0,1)
      
      #### model for marked individuals
      
      for (i in 1:mmax) #individual i
      {
      for(j in 1:J) #location j
      {
      D2[g,i,j] <- (S[g,i,1]-X[j,1])^2 + (S[g,i,2]-X[j,2])^2 #pythagorian theorem, distance to activity center
      lam[g,i,j] <- lam0*exp(-D2[g,i,j]/(2*sigma^2)) #lambda accounting for distance from activity center
      y[g,i,j]~dpois(lam.effm[g,i,j]*K[g,j])  #model accumulated counts across K for marked individuals and add for whether a camera is functioning or not
      lam.effm[g,i,j] <- lam[g,i,j]*zm[g,i] #lambda accounting for existence of individual
      }
      zm[g,i] ~ dbern(psim[g]) 
      S[g,i,1] ~ dunif(xlims[1], xlims[2]) #activity center x coord
      S[g,i,2] ~ dunif(ylims[1], ylims[2]) #activity center y coord

      }
      }
      for (g in 1:n.sess)
      {
      N[g] <- sum(zm[g,1:mmax])
      D[g] <- N[g]/S.Area
      }
    
      #### priors across grids
      lam0~dunif(0,5)
      sigma~dunif(0,100) #smoothing parameter; might need to adjust these based on what output looks like, I'm setting at approx 3*mmdm

      } #end model description
      
      ", fill = TRUE, file="SCR_multigrid.txt")
sink()

p <- function(sp.cr, n.chains, n.adapt, n.iter, n.burnin, n.cores)
{
  #load jags args
  load(paste("RData/CR-SCR/SCR_allgrids_jagsargs_", sp.cr, ".RData", sep=""))
  # n.chains=3; n.adapt=50; n.iter=250; n.burnin=50; n.cores=3 #test parameters
  n.chains=3; n.adapt=500; n.iter=2500; n.burnin=500; n.cores=3
  
  print(paste(n.chains, n.adapt, n.iter, n.burnin, n.cores))
  
  mod.SCR_B <- jags(data=data.SCR_B, inits=inits, parameters.to.save = params, model.file="Models/SCR_multigrid.txt", n.chains=n.chains,
                    n.adapt=n.adapt, n.iter=n.iter, n.burnin=n.burnin, n.thin=1, parallel=T, n.cores=n.cores)
  D.mode <- posterior.mode(as.mcmc(mod.SCR_B[[1]]$D))
  N.mode <- posterior.mode(as.mcmc(mod.SCR_B[[1]]$N))
  
  SCR_B <- list(sp.cr=sp.cr, mod=mod.SCR_B, Dhat=D.mode, Nhat=N.mode, SArea=data.SCR_B$S.Area)
  return(SCR_B)
}

pema <- p(sp.cr="PEMA", n.chains=3, n.adapt = 500, n.iter=2500, n.burnin=500, n.cores=3)
save(pema, paste(path.local, "RData/SCR_allgrids_", pema$sp.cr, ".RData", sep=""))

glsa <- p(sp.cr="GLSA", n.chains=3, n.adapt = 500, n.iter=2500, n.burnin=500, n.cores=3)
save(glsa, paste(path.local, "RData/SCR_allgrids_", glsa$sp.cr, ".RData", sep=""))

tato <- p(sp.cr="TATO", n.chains=3, n.adapt = 500, n.iter=2500, n.burnin=500, n.cores=3)
save(tato, paste(path.local, "RData/SCR_allgrids_", tato$sp.cr, ".RData", sep=""))

save(pema, glsa, tato, file="RData/SCR_allgrids_2023-03-10.RData")

############
#load jags output
load(paste(path.local, "RData/SCR_allgrids_TATO.RData", sep=""))

# output <- pema
# output <- glsa
# output <- tato
# 
# out <- capture.output(output)
# cat(out, file=paste("Output_SCR/", output$sp.cr, "_SCR-B_allgrids_", Sys.Date(), ".txt", sep=""), sep="\n", append=F)

# plot(output$mod)

require(MCMCvis)

s.all <- NULL
s2.all <- NULL
estimates <- function(model)
{
  s <- MCMCsummary(model$mod, params=c("lam0","sigma","psim","N","D"), round=3)
  s$param <- row.names(s)
  s$model <- "SCR_allgrids"
  s$sp <- model$sp.cr
  
  s.all <<- rbind(s.all, s)
  
  MCMCtrace(model$mod, params=c("lam0","sigma","psim","N"), type="density", filename=paste("Figures_SCR/SCR-B_density_allgrids_", model$sp.cr,".pdf", sep=""), open_pdf=F)
  MCMCtrace(model$mod, params=c("lam0","sigma","psim","N"), type="trace", filename=paste("Figures_SCR/SCR-B_trace_allgrids_", model$sp.cr,".pdf", sep=""), open_pdf=F)
  
  s2 <- data.frame(dhat=model$Dhat, nhat=model$Nhat, model="SCR_allgrids", sp=model$sp.cr, grid=1:length(model$Dhat))
  s2.all <<- rbind(s2.all, s2)
}

estimates(pema)
estimates(glsa)
estimates(tato)

write.table(s.all, file="Output/SCR-B_allgrids.txt", sep=",")
write.table(s2.all, file="Output/SCR-B_allgrids_mode.txt", sep=",")

# #plot mcmc chains
# jagsout.mcmc <- as.mcmc.list(output$mod$samples)
# p <- xyplot(jagsout.mcmc[,c("lam0", "sigma", "psim", "N")]) #plots traceplots in one figure
# dp <- densityplot(jagsout.mcmc[,c("lam0", "sigma", "psim", "N")]) #density plot of posterior
# 
# tiff(filename=paste("Figures_SCR/SCR-B_trace_", output$sp.cr, "_SM-0", g, "_", Sys.Date(), ".tiff", sep=""), height=6, width=6, res=300, units="in", compression="lzw")
# print(p)
# dev.off()
# 
# tiff(filename=paste("Figures_SCR/SCR-B_density_", output$sp.cr, "_SM-0", g, "_", Sys.Date(), ".tiff", sep=""), height=6, width=6, res=300, units="in", compression="lzw")
# print(dp)
# dev.off()
