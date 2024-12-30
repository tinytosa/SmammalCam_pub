
#multi grid, pooled data, capture-recapture spatial count models

##############################
#load packages
require(MCMCglmm)
require(MCMCvis)

############
#create objects to store all info from nimble output

s.all <- NULL
s2.all <- NULL
estimates <- function(model, priors)
{
  s <- MCMCsummary(model$mod, params=c("lam0","sigma","psi","N","D"), round=3)
  s$param <- row.names(s)
  s$model <- paste("CR-SC_", priors, "_multigrid", sep="")
  s$sp <- model$sp.cr
  s$g <- as.numeric(substr(s$param, start=nchar(s$param)-1, stop=nchar(s$param)-1))
  
  s.all <<- rbind(s.all, s)
  
  MCMCtrace(model$mod, params=c("lam0","sigma","psi","N","D"), type="density", ind=T, filename=paste("Figures_CR-SC/density_", priors,"_multigrid_", model$sp.cr,".pdf", sep=""))
  MCMCtrace(model$mod, params=c("lam0","sigma","psi","N","D"), type="trace", ind=T, filename=paste("Figures_CR-SC/trace_", priors, "_multigrid_", model$sp.cr,".pdf", sep=""), open_pdf=F)
  
  #modes of 
  modes <- data.frame(est=posterior.mode(model$mod))
  modes <- round(modes, digits=3)
  
  s2 <- data.frame(modes,  model=paste("CR-SC_", priors, "_multigrid", sep=""), sp=model$sp.cr)
  s2$param <- row.names(s2)
  s2.all <<- rbind(s2.all, s2)
}

############
#load jags/nimble output

load(file="RData/CR-SC_multigrid/CR-SC_multigrid_noinfo_PEMA.RData") #need to adjust grid numbers 6, 7, 8 -> 7, 8, 9
load(file="RData/CR-SC_multigrid/CR-SC_multigrid_noinfo_GLSA.RData")
load(file="RData/CR-SC_multigrid/CR-SC_multigrid_noinfo_TATO.RData")

estimates(pema, "noinfo")
#need to increase M for pema from M=500, M=1000 is good enough

estimates(tato, "noinfo")
#looks good

estimates(glsa, "noinfo")
#this is a mess! bad convergence but large enough M

rm(list = c('glsa','pema','tato'))

#############
#info
load(file="RData/CR-SC_multigrid/CR-SC_multigrid_info_PEMA.RData") #need to adjust grid numbers 6, 7, 8 -> 7, 8, 9
load(file="RData/CR-SC_multigrid/CR-SC_multigrid_info_GLSA.RData")
load(file="RData/CR-SC_multigrid/CR-SC_multigrid_info_TATO.RData")

estimates(pema, "info")
#need to increase M for pema from M=500, 1000, 2000 good enough (grid 7 looks like it's hitting the max)

estimates(tato, "info")
#need to increase M for tato from M=500, 1500 good enough, (grid 7, 8, 9 look like they're hitting the max)

estimates(glsa, "info")
#need to increase M for glsa from M=500, 1500 good enough (grid 6,7,8 look like they're hitting the max)

##########
#summarized in 4c_CR_spatialcount_separate_output
#don't need to do here

# write.table(s.all, file="Output_CR/CR-SC_multigrid_estimates.txt", sep=",", row.names=F)
# write.table(s2.all, file="Output_CR/CR-SC_multigrid_modes.txt", sep=",", row.names=F)
# 
# s <- merge(s.all, s2.all, by=c("param","model","sp"), all.x=T)

##########