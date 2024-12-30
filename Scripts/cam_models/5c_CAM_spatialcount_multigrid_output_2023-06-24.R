#multi grid, pooled data, camera trap spatial count models
#4 different t for each species: t=c(0, 15, 60, 1440)

##############################
#load packages
require(MCMCglmm)
require(MCMCvis)

############
#create objects to store all info from nimble output

estimates <- function(model, priors)
{
  s <- MCMCsummary(model$mod, params=c("lam0","sigma","psi","N","D"), round=3)
  s$param <- row.names(s)
  s$model <- paste("CAM-SC_", priors, "_multigrid", sep="")
  s$sp <- model$sp.cr
  s$g <- as.numeric(substr(s$param, start=nchar(s$param)-1, stop=nchar(s$param)-1))
  
  # s.all <<- rbind(s.all, s)
  
  MCMCtrace(model$mod, params=c("lam0","sigma","psi","N","D"), type="density", ind=T, filename=paste("Figures_CAM-SC/density_", priors,"_multigrid_", model$sp.cr,"_t", t, ".pdf", sep=""))
  MCMCtrace(model$mod, params=c("lam0","sigma","psi","N","D"), type="trace", ind=T, filename=paste("Figures_CAM-SC/trace_", priors, "_multigrid_", model$sp.cr,"_t", t, ".pdf", sep=""), open_pdf=F)
  
  #modes of 
  modes <- data.frame(est=posterior.mode(model$mod))
  modes <- round(modes, digits=3)
  
  s2 <- data.frame(modes,  model=paste("CAM-SC_", priors, "_multigrid", "_t", t, sep=""), sp=model$sp.cr)
  s2$param <- row.names(s2)
  s2$t <- t
  s2 <- merge(s2, s, by=c("sp","param"), all=T)
  s2.all <<- rbind(s2.all, s2)
}

############
#load jags/nimble output

# s.all <- NULL
s2.all <- NULL

for(t in c(0,15,60,1440))
{
  print(t)
  load(file = paste("RData/CAM-SC_multigrid_noinfo_PEMA_t",t,".RData", sep=""))
  estimates(pema, "noinfo")
}


for(t in c(0,15,60,1440))
{
  print(t)
  load(file = paste("RData/CAM-SC_multigrid_noinfo_TATO_t",t,".RData", sep=""))
  estimates(tato, "noinfo")
}


for(t in c(0,15,60,1440))
{
  print(t)
  load(file = paste("RData/CAM-SC_multigrid_noinfo_GLSA_t",t,".RData", sep=""))
  estimates(glsa, "noinfo")
}


###################
#info
for(t in c(0,15,60,1440))
{
  print(t)
  load(file = paste("RData/CAM-SC_multigrid_info_PEMA_t",t,".RData", sep=""))
  estimates(pema, "info")
}

#missing 15, 60
# for(t in c(0,15,60,1440))
for(t in c(0,1440))
{
  print(t)
  load(file = paste("RData/CAM-SC_multigrid_info_TATO_t",t,".RData", sep=""))
  estimates(tato, "info")
}

#missing 15, 60, 1440
# for(t in c(0,15,60,1440))
for(t in c(0))
{
  print(t)
  load(file = paste("RData/CAM-SC_multigrid_info_GLSA_t",t,".RData", sep=""))
  estimates(glsa, "info")
}

write.table(s2.all, file="Output_cam/output_CAM-SC_multigrid_estimates.txt", sep=",", row.names=T)
