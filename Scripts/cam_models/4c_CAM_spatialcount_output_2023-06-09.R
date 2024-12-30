###########################################
#spatial count models, info and no info, nimble MCMC run for each grid separately
#camera trap data
############################################

############
#load packages
require(MCMCglmm)
require(MCMCvis)

############
estimates <- function(model, priors, type = "separate")
{
  s <- MCMCsummary(model$mod, params=c("lam0","sigma","psi","N","D"), round=3)
  s$param <- row.names(s)
  s$model <- paste("CAM-SC_", priors, "_", type, sep="")
  if(type == "separate")
  {
    s$g <- model$g
  }
  if(type == "pool")
  {
    s$g <- "all"
    s[grep(s$param, pattern="\\["),]$g <- substr(s[grep(s$param, pattern="\\["),]$param, start=nchar(s[grep(s$param, pattern="\\["),]$param) - 1, stop=nchar(s[grep(s$param, pattern="\\["),]$param)-1)
  }
  s$sp <- model$sp.cr
  s$t <- model$t
  
  #print MCMC trace plots
  MCMCtrace(model$mod, params=c("lam0","sigma","psi","N","D"), type="density", ind=T, filename=paste("Figures_CAM-SC/density_", priors,"_", type, "_", model$sp.cr,"_", model$g, "_t", model$t, ".pdf", sep=""), open_pdf=F)
  MCMCtrace(model$mod, params=c("lam0","sigma","psi","N","D"), type="trace", ind=T, filename=paste("Figures_CAM-SC/trace_", priors, "_", type,"_", model$sp.cr,"_", model$g, "_t", model$t, ".pdf", sep=""), open_pdf=F)
  
  modes <- data.frame(mode=posterior.mode(as.mcmc.list(model$mod)))
  modes <- round(modes, digits=3)
  modes$param <- row.names(modes)
  
  s2.all <<- rbind(s2.all, merge(s, modes, by="param")) #combine mean and mode values
}

multigrid_estimates <- function(sp, priors="noinfo", type="separate")
{
  for(t in c(0,15,60,1440))
  {
    for(g in 1:8)
    {
      print(paste("t", t, "g", g))
      estimates(model=eval(parse(text=paste(sp, "_", g, "_t", t, sep=""))), priors=priors, type=type)
    }
  }
}



#################
#load RData
#no info
# rdata.list <- dir("RData/CAM-SC_separate", pattern="CAM-SC_separate_noinfo_", full.names=T)
rdata.list <- dir("RData/", pattern="CAM-SC_separate_noinfo_", full.names=T)

s2.all <- NULL

#updated on 10/1/2023
#pema grid 1-8, t=0,15,60,1440
load("RData/CAM-SC_separate/CAM-SC_separate_noinfo_PEMA_t0.RData") #t0, noinfo
load("RData/CAM-SC_separate/CAM-SC_separate_noinfo_PEMA_t1440.RData")  #t15, t60, t1440, noinfo
multigrid_estimates(sp="pema", priors="noinfo", type="separate")
rm(list=ls(pattern="pema"))

load("RData/CAM-SC_separate/CAM-SC_separate_info_PEMA_t1440.RData") #t0, t15, t60, t1440, info
multigrid_estimates(sp="pema", priors="info", type="separate")
rm(list=ls(pattern="pema"))


#tato
load("RData/CAM-SC_separate/CAM-SC_separate_noinfo_TATO_t0.RData") #t0
load("RData/CAM-SC_separate/CAM-SC_separate_noinfo_TATO_t1440.RData") #t15, t60, t1440
# sp <- "tato"
# priors <-"noinfo"
# type <- "separate"
# for(t in c("",15,60,1440))
# {
#   for(g in 1:8)
#   {
#     print(paste("t", t, "g", g))
#     estimates(model=eval(parse(text=paste(sp, "_", g, "_t", t, sep=""))), priors=priors, type=type)
#   }
# }
multigrid_estimates(sp="tato", priors="noinfo", type="separate")
rm(list=ls(pattern="tato"))

# load("RData/CAM-SC_separate/CAM-SC_separate_info_TATO_t1440.RData")
load("RData/CAM-SC_separate/CAM-SC_separate_info_TATO_t0.RData") #need to remove t15, t60, t1440 since those are all for noinfo, fine now after removing and saving
# rm(list=ls(pattern="_t15"))
# rm(list=ls(pattern="_t60"))
# rm(list=ls(pattern="_t1440"))
# save(list=ls(pattern="tato"), file="RData/CAM-SC_separate/CAM-SC_separate_info_TATO_t0.RData")
load("RData/CAM-SC_separate/CAM-SC_separate_info_TATO_t15.RData")
load("RData/CAM-SC_separate/CAM-SC_separate_info_TATO_t60.RData")
load("RData/CAM-SC_separate/CAM-SC_separate_info_TATO_t1440.RData")

multigrid_estimates(sp="tato", priors="info", type="separate")
rm(list=ls(pattern="tato"))

#glsa
load("RData/CAM-SC_separate/CAM-SC_separate_noinfo_GLSA_t1440.RData")
multigrid_estimates(sp="glsa", priors="noinfo", type="separate")
rm(list=ls(pattern="glsa"))

# load("RData/CAM-SC_separate/CAM-SC_separate_info_GLSA_SM-08_t1440.RData")
# load("RData/CAM-SC_separate/CAM-SC_separate_info_GLSA_t0.RData")
load("RData/CAM-SC_separate/CAM-SC_separate_info_GLSA_t1440.RData")
multigrid_estimates(sp="glsa", priors="info", type="separate")
rm(list=ls(pattern="glsa"))

################
write.table(s2.all, file="Output_cam/output_CAM-SC_separate_estimates.txt", sep=",", row.names = F)
