#spatial count model from capture recapture with small mammal data using informative priors
#each grid calculated separately

#adapted from D.Morin code 20190617_formatting example function
#built in figures for each animal's home range

#install.packages("MCMCglmm") #, repos='http://cran.us.r-project.org')
#install.packages("ggplot2", repos='http://cran.us.r-project.org', lib="/raid1/home/fw/tosam/R/powerpc64le-redhat-linux-gnu-library/3.5") #for ibm code

path.local <- "C:/Users/tosam/Documents/0_OSU/Dissertation/Data/SmallMammalCams/"

######
#install.packages("devtools")
#library(devtools)
#install_github("jaroyle/oSCR")

######
library(oSCR)
library(jagsUI) #have to install jags directly on your machine
library(MCMCglmm)
require(plyr) #for ddply
require(tidyr) #for separate
library(nimble)
library(coda)

###############
#load data
###
#capture recapture data
CR.data <- read.csv("WOOPS_2017_proofed_2019-06-10.csv", header=TRUE, sep=",") #all grids
CR.data$Loc <- paste(CR.data$Stn, CR.data$Stn.1, sep="")

CR.data <- CR.data[!CR.data$Loc=="NA",] #remove capture data rows with no trap location info
CR.data <- CR.data[!is.na(CR.data$SHER_DAY),] #remove data rows with no capture occasion number info

i <- ddply(CR.data, .(Species, Grid, Combined_Tag), nrow) #number of times each animal was captured
i <- ddply(ddply(CR.data, .(Species, Grid, Combined_Tag), nrow), .(Species, Grid, V1), nrow) #number of individuals captured on each grid

#write.table(table(i$Species, i$Grid), file="unique_inds.txt", sep=",")

####
#import trap data
traps.tmp <- read.csv("traplocs.txt", sep=",") #this includes both a large grid and a small grid for the CR data 

traps.all <- traps.tmp #if using all traps
traps.cam <- traps.tmp[traps.tmp$cam==1,] #remove half of small grid for where no cameras
traps.small <- traps.tmp[traps.tmp$Rlet > "H",] #removes traps from large grid, Sherman grid
traps.large <- traps.tmp[traps.tmp$Rlet < "I",] #removes traps from small grid, Tomahawk grid

#check each species and capture locations
check <- merge(CR.data, traps.all, by.x="Loc", by.y="RC", all.x=T)

table(check[check$Species %in% c("PEMA","TATO","GLSA"),c("ttype","Grid","Species")])

CR.data <- check[!(check$Species == "PEMA" & check$ttype == "tomahawk"),] #remove tomahawk traps for PEMA captures
CR.data <- CR.data[!(CR.data$Species != "PEMA" & CR.data$ttype == "sherman"),] #remove sherman traps for TATO and GLSA captures

table(CR.data[CR.data$Species %in% c("PEMA","TATO","GLSA"),c("ttype","Grid","Species")])

occ.trap <- data.frame(table(ddply(CR.data, .(Grid, SHER_DAY, ttype), nrow)[,c("Grid","ttype")])) #occasions per trap type per grid

#######
#formatting function
#######
format.fun<-function(data, traps, sp, grid, K, t=NULL)
{
  data$traps <- as.factor(paste(data$Stn,data$Stn.1, sep="")) #combine station letter and number
  #make edf, for indexing in package functions, need new trap numbers with no letters and in order match up correctly
  #length(levels(traps.tmp$RC))
  # traps$trap <- as.integer(traps$RC) #reassign trap numbers in CR data based on traps.tmp$traps
  traps$trap <- 1:nrow(traps) #assign trap number to traps
  
  lookup <- data.frame(trap.LET=traps$RC, trap.NUM=traps$trap) #create match key
  
  data$trap <- with(lookup,trap.NUM[match(data$traps,trap.LET)]) #match up to assigned numbers
  edf.data <- data.frame(session=as.integer(1), individual=data$Combined_Tag, occasion=data$SHER_DAY, TrapID=data$trap) #encounter data frame
  #edf.data <- data.frame(session=as.integer(1), individual=as.integer(data$Combined_Tag),occasion=1,TrapID=data$trap) #encounter data frame
  #I only have data for one session, but if choose to process all sessions for multi-strata model change up session, but you'll want it to be a vector of integers 1:session to work best 
  tdf.data <- list(data.frame(TrapID=traps[order(traps$trap),]$trap, Xcoord=traps[order(traps$trap),]$posx,Ycoord=traps[order(traps$trap),]$posy)) #trap data frame
  #,optional trap deployment matrix)) #replace with trap deployment file
  #if doing multi-strata, repeat the trap matrix above as the tdf list for length = number of strata
  scr.data <- data2oscr(edf=edf.data, sess.col=1, id.col=2, occ.col=3, trap.col=4, tdf=tdf.data, K = K, ntraps=c(nrow(traps)),
                        #scr.data <- data2oscr(edf=edf.data,sess.col=1,id.col=2,occ.col=3,trap.col=4,tdf=tdf.data,K = c(1),ntraps=c(nrow(traps)),
                        remove.extracaps = FALSE) #counts for Poisson
  sf <- scr.data$scrFrame
  sf
  tiff(filename=paste(path.local, "Figures_CR-SC/", sp,"_SM-", 0, grid, "_cr.tiff", sep=""), width=5, height=5, units="in", res=300, compression="lzw")
  print(plot(sf))
  dev.off()
  ssDF <- make.ssDF(scrFrame=sf, buffer=round(2*sf$mmdm,0), res=0.1*sf$mmdm, cont.cov = F, fact.cov = F) #used res=20 since trap res was 10 at fine scale
  plot(ssDF)
  points(traps$posx,traps$posy)
  out.formatted <- list(sf=sf,ssDF=ssDF,res=0.1*sf$mmdm)
  return(out.formatted)
}

##############

# g <- 1
# M = 500
# sp.cr <- "PEMA"
# traps <- traps.small
# 
# sp.cr <- "TATO"
# traps <- traps.large


p <- function(g, sp.cr, traps, M)
{
  print(paste(g, sp.cr))
  data.cr <- CR.data[which(CR.data$Species==sp.cr & CR.data$Grid==g),]
  data.cr <- data.cr[data.cr$Loc %in% traps$RC,]
  K <- occ.trap[occ.trap$Grid == g & occ.trap$ttype %in% unique(traps$ttype),]$Freq
  if(nrow(data.cr) < 2)
  {
    print(paste("not enough", sp.cr, "individuals/captures to format grid", g, Sys.Date()))
  }
  else
  {
    #####
    #set up unmarked model, can include non-ear tagged individuals
    data <- format.fun(data=data.cr, traps=traps, sp=sp.cr, grid=g, K=K) #switch to traps = large if the tomahawks, switch traps= small if sherman
  }
  
  yaug <- matrix(0, nrow=M, ncol=ncol(data$sf$caphist[[1]])) #formatting 0s for data augmentation
  yaug[1:nrow(data$sf$caphist[[1]]), ] <- apply(data$sf$caphist[[1]], 1:2, sum)  #and adding back in real individuals
  
  #n <- rowSums(yaug) #collapse individual encounter histories to number of recaps per individual
  n <- rowSums(colSums(data$sf$caphist[[1]]))
  S.Area <- nrow(data$ssDF[[1]])*data$res^2 #this in in meters^2
  
  K <- c(rep(occ.trap[occ.trap$Grid == g & occ.trap$ttype == "tomahawk",]$Freq, nrow(traps[traps$ttype == "tomahawk",])), 
         rep(occ.trap[occ.trap$Grid == g & occ.trap$ttype == "sherman",]$Freq, nrow(traps[traps$ttype == "sherman",])))
  
    ############
    #For JAGS
    # data.UNM <- list(M=M, J=ncol(data$sf$caphist[[1]]),
    #                  #K=rep(data$sf$occasions, times=ncol(data$sf$caphist[[1]])),
    #                  K=K,
    #                  n=n,
    #                  X=data$sf$traps[[1]],
    #                  xlims=c(min(data$ssDF[[1]]$X), max(data$ssDF[[1]]$X)),
    #                  ylims=c(min(data$ssDF[[1]]$Y), max(data$ssDF[[1]]$Y)),
    #                  S.Area=S.Area)
    # 
    # inits <- function(){list(lam0=runif(1), psi=runif(1), z=rep(1,data.UNM$M),
    #                          Su=cbind(data$ssDF[[1]]$X,data$ssDF[[1]]$Y)[sample(nrow(cbind(data$ssDF[[1]]$X,data$ssDF[[1]]$Y)),
    #                                                                  size=data.UNM$M,replace=TRUE),])} #sampling with replacement
    # params <- c("lam0", "sigma", "psi", "N", "D")
    # 
    # mod.UNM.jags <- jags(data=data.UNM, inits=inits, parameters.to.save = params, model.file=model.file, n.chains=n.chains,
    #                 n.adapt=n.adapt, n.iter=n.iter, n.burnin=n.burnin, n.thin=1,
    #                 parallel = TRUE, n.cores = n.cores)
    # D.mode <- posterior.mode(as.mcmc(mod.UNM.jags[[1]]$D)) #*10000 #/ha built into model instead
    # N.mode <- posterior.mode(as.mcmc(mod.UNM.jags[[1]]$N))
    # out.UNM<-list(mod=mod.UNM, Dhat=D.mode, Nhat=N.mode, SArea=S.Area)
    
    ############
    #For Nimble
    data.UNM <- list(n=n, X=data$sf$traps[[1]], K=K, 
                     xlims=c(min(data$ssDF[[1]]$X), max(data$ssDF[[1]]$X)),
                     ylims=c(min(data$ssDF[[1]]$Y), max(data$ssDF[[1]]$Y)))
    J <- ncol(data$sf$caphist[[1]]) #number of traps
    monitors <- c("lam0", "sigma", "psi", "N", "D") #same as params
    inits <- function(){list(lam0=runif(1), sigma=runif(1), psi=runif(1), z=rep(1, M),
                             #yu=matrix(rep(1, M*ncol(data$sf$caphist[[1]])), nrow=M, ncol=ncol(data$sf$caphist[[1]])),
                             Su=cbind(data$ssDF[[1]]$X, data$ssDF[[1]]$Y)[sample(nrow(cbind(data$ssDF[[1]]$X,data$ssDF[[1]]$Y)), size=constants$M, replace=T),])} #sampling with replacement
    
    save(data.UNM, M, J, S.Area, monitors, inits, sp.cr, file=paste("RData/CR-SC_separate_nimbleargs_", sp.cr,"_SM-0", g,".RData", sep=""))
}

# p(g=1, sp.cr="PEMA", traps=traps.small, M = 500)
# p(g=1, sp.cr="TATO", traps=traps.large, M = 500)

for(sp.cr in c("TATO","GLSA"))
{
  for(g in 1:9)
  {
    try(p(sp.cr=sp.cr, g=g, traps=traps.large, M=3000))
  }
}

for(sp.cr in c("PEMA"))
{
  for(g in c(1:5,7:9))
  {
    try(p(sp.cr=sp.cr, g=g, traps=traps.small, M=3000))
  }
}

#fix grid 8, TATO, K should be 6 not 7
sp.cr <- "TATO"
g <- 8
load(file=paste("RData/CR-SC_separate_nimbleargs_", sp.cr, "_SM-0", g,".RData", sep=""))
data.UNM$K <- rep(6, 64)
save(data.UNM, M, J, S.Area, monitors, inits, sp.cr, file=paste("RData/CR-SC_separate_nimbleargs_", sp.cr,"_SM-0", g,".RData", sep=""))
