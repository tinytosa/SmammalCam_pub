############################################
#spatial count models, info and no info, nimble MCMC run for each grid separately
############################################

############
require(MCMCglmm)
require(MCMCvis)

############
#load jags/nimble output
############

# model <- pema_1
# priors <- "info"

# s.all <- NULL
s2.all <- NULL

estimates <- function(model, priors, type = "separate")
{
  s <- MCMCsummary(model$mod, params=c("lam0","sigma","psi","N","D"), round=3)
  s$param <- row.names(s)
  s$model <- paste("CR-SC_", priors, "_", type, sep="")
  if(type == "separate")
  {
    s$g <- model$g
  }
  if(type == "pool")
  {
    s$g <- "all"
    s[grep(s$param, pattern="\\["),]$g <- substr(s[grep(s$param, pattern="\\["),]$param, start=nchar(s[grep(s$param, pattern="\\["),]$param) - 1, stop=nchar(s[grep(s$param, pattern="\\["),]$param)-1)
    
    s.area <- area[area$sp == model$sp.cr,] #fix for incorrect area
    s[grep(row.names(s), pattern="D"),]$mean <- s[grep(row.names(s), pattern="N"),]$mean/s.area$area*10000
    s[grep(row.names(s), pattern="D"),]$sd <- s[grep(row.names(s), pattern="N"),]$sd/s.area$area*10000
    s[grep(row.names(s), pattern="D"),]$`2.5%` <- s[grep(row.names(s), pattern="N"),]$`2.5%`/s.area$area*10000
    s[grep(row.names(s), pattern="D"),]$`50%` <- s[grep(row.names(s), pattern="N"),]$`50%`/s.area$area*10000
    s[grep(row.names(s), pattern="D"),]$`97.5%` <- s[grep(row.names(s), pattern="N"),]$`97.5%`/s.area$area*10000
  }
  s$sp <- model$sp.cr
  
  # s.all <<- rbind(s.all, s)
  
  #print MCMC trace plots
  MCMCtrace(model$mod, params=c("lam0","sigma","psi","N","D"), type="density", ind=T, filename=paste("Figures_CR-SC/density_", priors,"_", type, "_", model$sp.cr,"_", model$g, ".pdf", sep=""), open_pdf=T)
  # MCMCtrace(model$mod, params=c("lam0","sigma","psi","N","D"), type="trace", ind=T, filename=paste("Figures_CR-SC/trace_", priors, "_", type,"_", model$sp.cr,"_", model$g,".pdf", sep=""), open_pdf=F)
  
  modes <- data.frame(mode=posterior.mode(model$mod))
  modes <- round(modes, digits=3)
  modes$param <- row.names(modes)
  
  if(type == "pool")
  {
    modes[grep(row.names(modes), pattern="D"),]$mode <- modes[grep(row.names(modes), pattern="N"),]$mode/s.area$area*10000
  }
  
  s2.all <<- rbind(s2.all, merge(s, modes, by="param")) #combine mean and mode values
  # s2 <- data.frame(dhat=modes["D",], nhat=modes["N",], model=paste("CR-SC_", priors, "_", type, sep=""), sp=model$sp.cr, grid=model$g)
  # s2.all <<- rbind(s2.all, s2)
}

#################
#load RData
#no info
load("RData/CR-SC_separate/CR-SC_separate_noinfo_PEMA.RData")
estimates(pema_1, priors="noinfo")
estimates(pema_2, priors="noinfo")
estimates(pema_3, priors="noinfo")
estimates(pema_4, priors="noinfo")
estimates(pema_5, priors="noinfo")
estimates(pema_7, priors="noinfo")
estimates(pema_8, priors="noinfo")
estimates(pema_9, priors="noinfo")

#M=500 is not enough for PEMA noinfo
#M=1000 is enough, but no convergence

#
load("RData/CR-SC_separate/CR-SC_separate_noinfo_GLSA.RData") #M=500
estimates(glsa_1, priors="noinfo")
estimates(glsa_2, priors="noinfo")
estimates(glsa_3, priors="noinfo")
estimates(glsa_4, priors="noinfo")
estimates(glsa_5, priors="noinfo")
estimates(glsa_6, priors="noinfo")
estimates(glsa_7, priors="noinfo")
estimates(glsa_8, priors="noinfo")
estimates(glsa_9, priors="noinfo")

#M=500 is enough

#
load("RData/CR-SC_separate/CR-SC_separate_noinfo_TATO.RData") #M=500
estimates(tato_1, priors="noinfo")
estimates(tato_2, priors="noinfo")
estimates(tato_3, priors="noinfo")
estimates(tato_4, priors="noinfo")
estimates(tato_5, priors="noinfo")
estimates(tato_6, priors="noinfo")
estimates(tato_7, priors="noinfo")
estimates(tato_8, priors="noinfo")
estimates(tato_9, priors="noinfo")
#SArea = 334374.4 m^2

#M=500 is enough

rm(list=ls(pattern="pema"))
rm(list=ls(pattern="glsa"))
rm(list=ls(pattern="tato"))

###################################
#load RData
#info
load("RData/CR-SC_separate/CR-SC_separate_info_PEMA.RData")
estimates(pema_1, priors="info")
estimates(pema_2, priors="info")
estimates(pema_3, priors="info")
estimates(pema_4, priors="info")
estimates(pema_5, priors="info")
estimates(pema_7, priors="info")
estimates(pema_8, priors="info")
estimates(pema_9, priors="info")
#M=500, 1000, 2000 is not enough for PEMA info

load("RData/CR-SC_separate/CR-SC_separate_info_GLSA.RData")
estimates(glsa_1, priors="info")
estimates(glsa_2, priors="info")
estimates(glsa_3, priors="info")
estimates(glsa_4, priors="info")
estimates(glsa_5, priors="info")
estimates(glsa_6, priors="info")
estimates(glsa_7, priors="info")
estimates(glsa_8, priors="info")
estimates(glsa_9, priors="info")
#M=500, 1000, 2000 is not enough for GLSA info

load("RData/CR-SC_separate/CR-SC_separate_info_TATO.RData")
estimates(tato_1, priors="info")
estimates(tato_2, priors="info")
estimates(tato_3, priors="info")
estimates(tato_4, priors="info")
estimates(tato_5, priors="info")
estimates(tato_6, priors="info")
estimates(tato_7, priors="info")
estimates(tato_8, priors="info")
estimates(tato_9, priors="info")
#M=500, 1000, 2000 is not enough for TATO info

rm(list=ls(pattern="pema"))
rm(list=ls(pattern="glsa"))
rm(list=ls(pattern="tato"))

###################################
#load RData
#need to fix density estimates since SArea was wrong
area <- read.table("Output_CR/multigrid_area.txt", sep=",", header=T)

#no info, pooled (multigrid)
load("RData/CR-SC_multigrid/CR-SC_multigrid_noinfo_PEMA.RData")
estimates(pema, priors="noinfo", type="pool") #not converged, but M=500, 1000 is enough

load("RData/CR-SC_multigrid/CR-SC_multigrid_noinfo_TATO.RData")
estimates(tato, priors="noinfo", type="pool") #good with M=500

load("RData/CR-SC_multigrid/CR-SC_multigrid_noinfo_GLSA.RData")
estimates(glsa, priors="noinfo", type="pool") #not converged, but M=500 is enough

rm(list=ls(pattern="pema"))
rm(list=ls(pattern="glsa"))
rm(list=ls(pattern="tato"))

###################################
#load RData
#info, pooled (multigrid)

load("RData/CR-SC_multigrid/CR-SC_multigrid_info_PEMA.RData")
estimates(pema, priors="info", type="pool") #not converged, but M=2000 is good enough for grid 2,7 (actually 8)

load("RData/CR-SC_multigrid/CR-SC_multigrid_info_TATO.RData")
estimates(tato, priors="info", type="pool") #not converged, but M=500 is not enough for grid 2,3,5,7 (actually 8),8 (actually 9)
#running M=1500 6/9/2023, same as M=500

load("RData/CR-SC_multigrid/CR-SC_multigrid_info_GLSA.RData")
estimates(glsa, priors="info", type="pool") #not converged, but M=1000 is not enough for grid 2,7 (actually 8)
#running M=1500 on 6/9/2023, same as M=500

###################################
# write.table(s.all, file="Output_CR/CR-SC_separate_estimates.txt", sep=",", row.names=F) #don't need this anymore, added mean and mode to s2.all
write.table(s2.all, file="Output_CR/CR-SC_estimates.txt", sep=",", row.names=F)

# fix grid numbers for 7,8,9 for PEMA
s2.all <- read.table("Output_CR/CR-SC_estimates.txt", sep=",", header=T)

# s2.all <- s2.all %>% mutate_if(is.numeric, round, digits=2) %>% mutate_if(is.numeric, as.character) #round to 2 decimal points, this cuts off 0s at the end of the number
s2.all <- s2.all %>% mutate_if(is.numeric, sprintf, fmt='%#.2f') #round to 2 decimal points
s2.all <- s2.all[grep(s2.all$param, pattern="D"),] #only keep density estimates
s2.all[s2.all$sp == "PEMA",]$g <- rep(c(1:5,7:9), 4) #correct grid numbers for pema
s2.all$mode_range <- paste(s2.all$mode, " (", s2.all$X2.5.," - ", s2.all$X97.5., ")", sep="") #convert to mode (range)
s2.all <- s2.all[s2.all$Rhat < 1.1,] #only keep values for models that converged

s2.all$sp <- factor(s2.all$sp, levels=c("PEMA","TATO","GLSA"), labels=c("mouse","chipmunk","flying squirrel"))

s2.all <- s2.all[order(s2.all$sp, s2.all$g),]

#convert to wide
s2.wide <- pivot_wider(s2.all[,c("model","g","sp","mode_range")], names_from=model, values_from=mode_range)
s2.wide <- s2.wide[,c('sp','g','CR-SC_noinfo_separate','CR-SC_info_separate','CR-SC_noinfo_pool','CR-SC_info_pool')]
write.csv(s2.wide, file="Output_CR/all_CR-SC.csv")

##########################
# results <- as.mcmc.list(mod.UNM)
# modes <- posterior.mode(results)
# rhat <- gelman.diag(results, multivariate = F)
# # posterior.mode(results$chain1)
# # posterior.mode(results$chain2)
# # posterior.mode(results$chain3)
# out.UNM <- list(mod=mod.UNM, results=results, modes=modes, rhat=rhat)
# 
# out <- capture.output(summary(UNM_B.cr$mod), UNM_B.cr$modes, UNM_B.cr$rhat, UNM_B.cr$SArea)
# cat(out, file=paste("Output_UNM-cr/nimble_", sp.cr, "_SM-0", g, "_UNM-B.cr_", Sys.Date(), ".txt", sep=""), sep="\n", append=F)
# 
# tp <- xyplot(UNM_B.cr$results) #plots traceplots in one figure
# dp <- densityplot(UNM_B.cr$results) #density plot of posterior
# 
# tiff(filename=paste("Figures_UNM-cr/outnimble_trace_", sp.cr, "_SM-0", g, "_", Sys.Date(),".tiff", sep=""), height=6, width=6, res=300, units="in", compression="lzw")
# print(tp)
# dev.off()
# 
# tiff(filename=paste("Figures_UNM-cr/outnimble_densityplot_", sp.cr, "_SM-0", g, "_", Sys.Date(), ".tiff", sep=""), height=6, width=6, res=300, units="in", compression="lzw")
# print(dp)
# dev.off()
