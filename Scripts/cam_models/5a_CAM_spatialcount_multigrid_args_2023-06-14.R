

######
#load packages

library(oSCR)
library(jagsUI) #have to install jags directly on your machine
library(MCMCglmm)
require(plyr) #for ddply
require(tidyr) #for separate
library(nimble)
library(coda)

######
#load data
#camera data
cam.data <- read.table("SMcamera_tags_daynum_2020-05-04.txt", sep=",") #only inlcudes data for deer mice, flying squirrels, townsends chipmunk, and voles
cam.data <- cam.data[cam.data$num %in% c(1:8),]
cam.data <- cam.data[!cam.data$grid == 0,]

#import trap data
traps.tmp <- read.csv("traplocs.txt", sep=",") #this includes both a large grid and a small grid for the CR data 

#as.integer for characters stopped working
# traps.small <- traps.tmp[as.integer(traps.tmp$Rlet) > 8,] #removes traps from large grid
# traps.large <- traps.tmp[as.integer(traps.tmp$Rlet) < 9,] #removes traps from small grid
traps.small <- traps.tmp[traps.tmp$Rlet > "H",]
traps.large <- traps.tmp[traps.tmp$Rlet < "I",]

traps.all <- traps.tmp #if using all traps
traps.cam <- traps.tmp[traps.tmp$cam == 1,] #remove half of small grid for where no cameras

traps.both <- merge(traps.large, traps.small, by=c("posx","posy"), suffixes=c(".t",".s")) #find RC of tomahawk traps that correspond to sherman grid
traps.cam.small <- traps.cam[traps.cam$RC %in% traps.both$RC.t | traps.cam$RC %in% traps.small$RC,] #nested cameras along sherman grid
traps.cam.large <- traps.cam[traps.cam$ttype == "tomahawk",]

occ.trap <- data.frame(table(ddply(cam.data, .(grid, num), nrow)[,c("grid")]))
names(occ.trap) <- c("Grid","Freq")

traps.cam.all <- data.frame(RC=rep(traps.cam$RC, 8), grid=rep(1:8, each=80))
traps.cam.small <- traps.cam[traps.cam$RC %in% traps.both$RC.t | traps.cam$RC %in% traps.small$RC,] #nested cameras along sherman grid
traps.cam.large <- traps.cam[traps.cam$ttype == "tomahawk",]

plot(traps.cam.small$posx, traps.cam.small$posy)
plot(traps.cam.large$posx, traps.cam.large$posy)

#######
#formatting function
#######
format.fun <- function(data, traps, sp, t=NULL)
{
    #consolidate camera photos if they are t mins apart
    data$DT <- as.POSIXct(data$DateTimeOriginal, format="%Y-%m-%d %H:%M:%S", tz="GMT")
    data$diff <- c(0,round((data$DT[-1] - data$DT[-nrow(data)])/60, digits=0)) #time differences in minutes
    
    #function to number rows to determine which to consolidate
    n <- 1
    data$consolidate <- sapply(data$diff, function(x) {if(x > t ){n <<- n+1} 
      n})
    
    #consolidate if keywords and consolidate numbers are the same
    data <- ddply(data, .(Directory, grid, Loc, Rlet, Cnum, Keywords, num, consolidate), summarise, numphoto=length(DT), start=min(DT))
    
    # print(table(data[,c("Rlet","Cnum")]))
    # write.table(paste(Sys.Date(), "g=",g, "t=",t, sp.cr), file=paste(path.local, "consolidated_data.txt", sep=""), sep=",", append=T, col.names=F)
    # write.table(table(data[,c("Rlet","Cnum")]), file=paste(path.local, "consolidated_data.txt", sep=""), sep=",", append=T)
    
    data$traps <- as.factor(data$Loc)
    # traps$trap <- as.integer(as.factor(traps$RC)) #reassign trap numbers in CR data based on traps.tmp$traps
    traps$trap <- 1:nrow(traps) #assign trap number to traps
    
    lookup <- data.frame(trap.LET=traps$RC, trap.NUM=traps$trap) #create match key
    data$trap <- with(lookup, trap.NUM[match(data$traps, trap.LET)]) #match up to assigned numbers
    data$Combined_Tag <- 1
    
    edf.data <- data.frame(session=data$grid, individual=as.integer(data$Combined_Tag), #trying to fake oSCR out by making all one individual
                           occasion=as.integer(1), #using Poisson so can put all in same occasion, but if add in trap operation need to make full matrix
                           #occasion=as.integer(data$num), #separate detections by occasion
                           TrapID=data$trap)
    
    tdf.data <- data.frame(TrapID=traps[order(traps$trap),]$trap, Xcoord=traps[order(traps$trap),]$posx,Ycoord=traps[order(traps$trap),]$posy)
    unm.data <- data2oscr(edf=edf.data, 
                          sess.col=1, id.col=2, occ.col=3, trap.col=4,
                          tdf=rep(list(tdf.data), 8), #one per session 
                          K = rep(1,8), #one per session
                          ntraps=rep(nrow(traps),8), #one per session
                          remove.extracaps = FALSE) #counts for Poisson
    sf <- unm.data$scrFrame
    sf
    
    ssDF <- make.ssDF(scrFrame=sf, buffer=round(1.5*(max(sf$traps[[1]][,1])-min(sf$traps[[1]][,1])),0), #don't have mmdm so making 1.5*grid width
                      res=10, cont.cov = F, fact.cov = F) #used res=10m since trap res was 10 at fine scale
    #this will take a lot longer to run
    plot(ssDF)
    points(traps$posx,traps$posy)
    
    # K = matrix(data=rep(data$sf$occasions, times=ncol(data$sf$caphist[[1]])), nrow=length(data$sf$caphist), ncol=ncol(data$sf$caphist[[1]])) #matrix of [grid, number of occasions]
    
    K.all <- ddply(data, .(grid, num, Loc), nrow)
    K.all$occ <- 1
    traps.grids <- data.frame(rbind(traps, traps, traps, traps, traps, traps, traps, traps), grid=rep(1:8, each=nrow(traps)))
    K.all <- merge(K.all, traps.grids, by.x=c("Loc","grid"), by.y=c("RC","grid"), all.y=T) #make sure all trap locations are represented
    if(nrow(K.all[is.na(K.all$occ),]) > 0)
    {
      K.all[is.na(K.all$occ),]$occ <- 0
    }
    K.all <- ddply(K.all, .(grid, Loc, occ), nrow) #consolidate, creates one row for every trap location in each grid
    K.all <- K.all[,c("grid","Loc","occ")]
    K.all$occ <- K.all$occ * 8 #get number of days each camera was operational
    
    
    K.all$LocNum <- with(lookup, trap.NUM[match(K.all$Loc, trap.LET)]) #match up to assigned numbers
    K.all <- K.all[order(K.all$grid, K.all$LocNum),] #make order of K.all match n below
    K <- matrix(K.all$occ, nrow=8, ncol=nrow(traps), byrow=T)
    
    out.formatted <- list(sf=sf, ssDF=ssDF, res=10, K=K)
    
  return(out.formatted)
}

p <- function(traps, sp.cam, sp.cr, M)
{
  # traps <- traps.cam.small
  # sp.cam <- "DEER MOUSE"
  # sp.cr <- "PEMA"
  # M <- 10
  # t <- 1440
  
  print(paste(sp.cam, "M", M))
  data.cam <- cam.data[which(cam.data$Keywords==sp.cam),] #subset data for species and grid
  data.cam <- data.cam[data.cam$Loc %in% traps$RC,] #subset so that only data from cameras in "traps"
  
  #######
  for(t in c(0,15,60,1440))
  {
    print(t)
    cam <- format.fun(data=data.cam, traps=traps, sp=sp.cam, t=t)
    
    n <- t(sapply(cam$sf$caphist, function(x) {rowSums(colSums(x))})) # matrix of [grid, number of captures per trap loc] summed over all capture occasions
    print(n)
    
    S.Area <- nrow(cam$ssDF[[1]])*cam$res^2 #this in in meters
    print(S.Area)
    
    K <- cam$K
    print(K)
    
    J <- ncol(cam$sf$caphist[[1]]) #number of traps
    
    ############
    #For Nimble
    ############
    locs <- cbind(cam$ssDF[[1]]$X, cam$ssDF[[1]]$Y)
    data.UNM <- list(n=n, X=cam$sf$traps[[1]], 
                     xlims=c(min(cam$ssDF[[1]]$X), max(cam$ssDF[[1]]$X)),
                     ylims=c(min(cam$ssDF[[1]]$Y), max(cam$ssDF[[1]]$Y)))
    monitors <- c("lam0", "sigma", "psi", "N", "D") #same as params
    inits <- function()
      {
      list(lam0=runif(1), sigma=runif(1,0,100),
           psi=runif(length(cam$sf$caphist)), #psi for each grid (8)
           z=matrix(data=rbinom(M,1,0.5), nrow=length(cam$sf$caphist), ncol=M), #pull from binomial distribution, 8 rows x M col
           #yu=matrix(rep(1, M*ncol(data$sf$caphist[[1]])), nrow=M, ncol=ncol(data$sf$caphist[[1]])),
           Su=array(data=locs[sample(nrow(locs), size=M*length(cam$sf$caphist), replace=T),], dim=c(length(cam$sf$caphist), M, 2)))
      } #dim = num grids x M rows x 2 for location of each M
    # save(cam, data.UNM, M, K, J, t, S.Area, monitors, inits, sp.cam, sp.cr, file=paste("RData/CAM-SC_multigrid_nimbleargs_", sp.cr,"_test.RData", sep=""))
    save(cam, data.UNM, M, K, J, t, S.Area, monitors, inits, sp.cam, sp.cr, file=paste("RData/CAM-SC_multigrid_nimbleargs_", sp.cr,"_t",t,".RData", sep=""))
  }
}

p(sp.cam="DEER MOUSE", sp.cr="PEMA", traps=traps.cam.small, M=2000)
p(sp.cam="TOWNSENDS CHIPMUNK", sp.cr="TATO", traps=traps.cam.large, M=2000)
p(sp.cam="FLYING SQUIRREL", sp.cr="GLSA", traps=traps.cam.large, M=2000)