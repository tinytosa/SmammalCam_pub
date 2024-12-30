
#######################
#load packages
library(oSCR)
library(jagsUI) #have to install jags directly on your machine
library(MCMCglmm)
require(plyr) #for ddply
require(tidyr) #for separate
library(nimble)
library(coda)

#######################
#laod data
path.local <- "C:/Users/tosam/Documents/0_OSU/Dissertation/Data/SmallMammalCams/"

CR.data<-read.csv("WOOPS_2017_proofed_2019-06-10.csv", header=TRUE, sep=",") #all grids
CR.data$Loc <- paste(CR.data$Stn, CR.data$Stn.1, sep="")

CR.data <- CR.data[!CR.data$Loc=="NA",] #remove data rows with no trap location info
CR.data <- CR.data[!is.na(CR.data$SHER_DAY),] #remove data rows with no capture occasion number info

i <- ddply(CR.data, .(Species, Grid, Combined_Tag), nrow)
i <- ddply(ddply(CR.data, .(Species, Grid, Combined_Tag), nrow), .(Species, Grid, V1), nrow)

#import trap data
traps.tmp<-read.csv("traplocs.txt", sep=",") #this includes both a large grid and a small grid for the CR data 

traps.all <- traps.tmp #if using all traps
traps.cam <- traps.tmp[traps.tmp$cam==1,] #remove half of small grid for where no cameras
traps.small <- traps.tmp[traps.tmp$Rlet > "H",] #removes traps from large grid, Sherman grid
traps.large <- traps.tmp[traps.tmp$Rlet < "I",] #removes traps from small grid, Tomahawk grid

check <- merge(CR.data, traps.all, by.x="Loc", by.y="RC", all.x=T)

table(check[check$Species %in% c("PEMA","TATO","GLSA"),c("ttype","Grid","Species")])

CR.data <- check[!(check$Species == "PEMA" & check$ttype == "tomahawk"),] #remove tomahawk traps for PEMA captures
CR.data <- CR.data[!(CR.data$Species != "PEMA" & CR.data$ttype == "sherman"),] #remove sherman traps for TATO and GLSA captures

table(CR.data[CR.data$Species %in% c("PEMA","TATO","GLSA"),c("ttype","Grid","Species")])

CR.data <- check[!(check$Species == "PEMA" & check$ttype == "tomahawk"),] #remove tomahawk traps for PEMA captures
CR.data <- CR.data[!(CR.data$Species != "PEMA" & CR.data$ttype == "sherman"),] #remove sherman traps for TATO and GLSA captures

table(CR.data[CR.data$Species %in% c("PEMA","TATO","GLSA"),c("ttype","Grid","Species")])


occ.trap <- data.frame(table(ddply(CR.data, .(Grid, SHER_DAY, ttype), nrow)[,c("Grid","ttype")])) #occasions per trap type per grid

#######
#formatting function
#######
format.fun<-function(data, traps, sp)
{ 
  data$traps <- as.factor(paste(data$Stn,data$Stn.1, sep="")) #combine station letter and number
  #make edf, for indexing in package functions, need new trap numbers with no letters and in order match up correctly
  #length(levels(traps.tmp$RC))
  traps$trap <- 1:nrow(traps) #assign trap number to traps
  
  lookup <- data.frame(trap.LET=traps$RC,trap.NUM=traps$trap) #create match key
  
  data$trap <- with(lookup,trap.NUM[match(data$traps,trap.LET)]) #match up to assigned numbers
  data <- data[order(data$Grid),] #sort data by grid
  
  edf.data <- data.frame(session=data$Grid, individual=data$Combined_Tag, occasion=data$SHER_DAY, TrapID=data$trap) #encounter data frame
  print(table(edf.data[,c("session","occasion")]))
  #I only have data for one session, but if choose to process all sessions for multi-strata model change up session, but you'll want it to be a vector of integers 1:session to work best 
  tdf.data <- data.frame(TrapID=traps[order(traps$trap),]$trap, Xcoord=traps[order(traps$trap),]$posx,Ycoord=traps[order(traps$trap),]$posy) #trap data frame
  #,optional trap deployment matrix)) #replace with trap deployment file
  #if doing multi-strata, repeat the trap matrix above as the tdf list for length = number of strata
  
  table(edf.data[,c("session","occasion")])
  if(sp == "PEMA") #since no grid 6 data
  {
    #use different values for session
    data$gnew <- data$Grid
    data[data$gnew == 7,]$gnew <- 6
    data[data$gnew == 8,]$gnew <- 7
    data[data$gnew == 9,]$gnew <- 8
    edf.data <- data.frame(session=data$gnew, individual=data$Combined_Tag, occasion=data$SHER_DAY, TrapID=data$trap) #encounter data frame
    
    edf.data <- edf.data[!(edf.data$session == 7 & edf.data$occasion == 4),] #only 1 capture in grid 8, session 4; weather?
    #no grid 6 data
    scr.data <- data2oscr(edf=edf.data,
                          sess.col=1, id.col=2, occ.col=3, trap.col=4, 
                          tdf=rep(list(tdf.data), 8), #one for each session
                          K = c(4,4,4,4,4,4,3,2), #one per session
                          ntraps=rep(nrow(traps), 8), #one per session
                          remove.extracaps = FALSE) #counts for Poisson
  }
  else
  {
    K.info <- unique(edf.data[,c("session","occasion")])
    # K.info <- K.info[order(K.info$occasion),]
    # K.info <- K.info[order(K.info$session),]
    K.info <- ddply(K.info, .(session), summarize, max.occ=max(occasion))
    scr.data <- data2oscr(edf=edf.data, 
                          sess.col=1, id.col=2, occ.col=3, trap.col=4, 
                          tdf=rep(list(tdf.data), 9),
                          #tdf=list(tdf.data, tdf.data, tdf.data, tdf.data, tdf.data, tdf.data, tdf.data, tdf.data, tdf.data), #one for each session
                          # K = c(8,8,8,4,8,6,8,7,4), #one per session
                          K = K.info$max.occ,
                          ntraps=rep(nrow(traps), 9), #one per session
                          remove.extracaps = FALSE) #counts for Poisson
  }
  
  sf <- scr.data$scrFrame
  sf
  # tiff(filename=paste(path.local, "Figures_SCR/", sp, "_cr.tiff", sep=""), width=5, height=5, units="in", res=300, compression="lzw")
  # print(plot(sf))
  # dev.off()
  ssDF <- make.ssDF(scrFrame=sf, buffer=round(3*sf$mmdm,0), res=0.1*sf$mmdm, cont.cov = F, fact.cov = F) #used res=20 since trap res was 10 at fine scale
  plot(ssDF)
  points(traps$posx,traps$posy)
  out.formatted <- list(sf=sf,ssDF=ssDF,res=0.1*sf$mmdm)
  return(out.formatted)
}

p <- function(sp.cr, traps, M)
{
  data.cr <- CR.data[which(CR.data$Species==sp.cr),]
  data.cr <- data.cr[data.cr$Loc %in% traps$RC,]
  
  if(nrow(data.cr) < 2)
  {
    print(paste("not enough", sp.cr, "individuals/captures to format grid", g, Sys.Date()))
  }
  else
  {
    #####
    #set up unmarked model, can include non-ear tagged individuals
    data <- format.fun(data=data.cr, traps=traps, sp=sp.cr) #switch to traps = large if the tomahawks, switch traps= small if sherman
  }
  
  #set up data
  y3d <- data$sf$caphist #extracting capture history array from oSCR scr frame
  S.Area <- nrow(data$ssDF[[1]])*data$res^2 #this in in meters
  print(S.Area)
  
  #create 3D matrix: grid x ninds x traps
  yaug <- array(data=0, dim=c(length(y3d), M, ncol(data$sf$caphist[[1]]))) #dim(num sessions, mmax, num traps)
  for(g in 1:length(y3d))
  {
    yaug[g, 1:nrow(y3d[[g]]),] <- apply(y3d[[g]], 1:2, sum) #check this
  }
  
  locs <- cbind(data$ssDF[[1]]$X, data$ssDF[[1]]$Y)
  n <- t(sapply(data$sf$caphist, function(x) {rowSums(colSums(x))})) # matrix of [grid, number of captures per trap loc] summed over all capture occasions
  K = matrix(data=rep(data$sf$occasions, times=ncol(data$sf$caphist[[1]])), nrow=length(data$sf$caphist), ncol=ncol(data$sf$caphist[[1]])) #matrix of [grid, number of occasions]
  print(data$sf$occasions)
  J = ncol(data$sf$caphist[[1]]) #number of traps
  
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
  #change this for multigrid
  ############
  data.UNM <- list(n=n, X=data$sf$traps[[1]], K=K,
                   xlims=c(min(data$ssDF[[1]]$X), max(data$ssDF[[1]]$X)),
                   ylims=c(min(data$ssDF[[1]]$Y), max(data$ssDF[[1]]$Y)))
  monitors <- c("lam0", "sigma", "psi", "N", "D") #same as params
  inits <- function(){list(lam0=runif(1), sigma=runif(1, 0, 200), psi=runif(length(data$sf$caphist), 0.1, 0.9),
                           # z=matrix(data=1, nrow=length(y3d), ncol=M), #for single grid, z=rep(1, M),
                           z=matrix(data=rbinom(1,1,0.5), nrow=length(y3d), ncol=M),
                           Su=array(data=locs[sample(nrow(locs), size=M*length(y3d), replace=T),], dim=c(length(y3d), M, 2))
  )} #sampling with replacement
  
  save(data.UNM, M, J, S.Area, monitors, inits, sp.cr, file=paste("RData/CR-SC_multigrid_nimbleargs_", sp.cr,".RData", sep=""))
}


p(M=2000, sp.cr="PEMA", traps=traps.small) #changed M from 500 to 1000
p(M=1500, sp.cr="TATO", traps=traps.large)
p(M=1500, sp.cr="GLSA", traps=traps.large)

area <- data.frame(sp=c("PEMA","TATO","GLSA"), area=c(55820.32, 567639, 565299.7))
write.table(area, file="Output_CR/multigrid_area.txt", sep=",", row.names=F)
