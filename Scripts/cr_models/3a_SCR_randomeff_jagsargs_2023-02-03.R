
#spatial capture recapture with small mammal data
#adapted from D.Morin code 20190617_formatting example function
#built in figures for each animal's home range
#hierarchical model with grid as a random effect

############################
#install.packages("MCMCglmm") #, repos='http://cran.us.r-project.org')
#install.packages("ggplot2", repos='http://cran.us.r-project.org', lib="/raid1/home/fw/tosam/R/powerpc64le-redhat-linux-gnu-library/3.5") #for ibm code

# path.local <- "C:/Users/tosam/Documents/0_OSU/Dissertation/Data/SmallMammalCams/"
path.local <- "C:/Users/tosam/Documents/0_OSU/Dissertation/Manuscripts/SmallMammalCameras/SmammalCam_pub/"
#path.local <- "/raid1/home/fw/tosam/SmallMammalCams/"
setwd(path.local)

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

#20190620 Notes from MARIE on grids:

# Columns Rlet and Cnum in the SMcamera file represent trap stations.
# Columns Stn and Stn in the WOOPS live trap file represent trap stations.
# Rlet and Stn letters belong on the y axis. Letters A-H should correspond to positions: 280,240,200,160,120,80,40,0 OR  rev(seq(from=0, by=40, length.out = 8)); Letters P-Y should correspond to positions: 160,150,140,130,120,110,100,90,80,70 OR seq(from=160, by=-10, length.out=10)
# Cnum and stn numbers belong on the x axis. Numbers for stations A-H 1-8 correspond to positions: 0,40,80,120,160,200,240,280 OR seq(0, 280, by=40). Numbers for stations P-Y 1-10 correspond to positions: 160,170,180,190,200,210,220,230,240,250 OR seq(160, by=10, length.out=10)

#load and see what each of these files look like.
CR.data<-read.csv(paste(path.local,"Data_Raw/WOOPS_2017_proofed_2019-06-10.csv", sep=""), header=TRUE, sep=",")
CR.data$Loc <- paste(CR.data$Stn, CR.data$Stn.1, sep="")
CR.data <- CR.data[!CR.data$Loc=="NA",] #remove data rows with no trap location info
CR.data <- CR.data[!is.na(CR.data$SHER_DAY),] #remove data rows with no capture occasion number info

i <- ddply(CR.data, .(Species, Grid, Combined_Tag), nrow)
i <- ddply(ddply(CR.data, .(Species, Grid, Combined_Tag), nrow), .(Species, Grid, V1), nrow)
#write.table(table(i$Species, i$Grid), file="unique_inds.txt", sep=",")

#import trap data
traps.tmp <- read.csv(paste(path.local, "Data_Raw/traplocs.txt", sep=""), sep=",") #this includes both a large grid and a small grid for the CR data 
#(different trap types catch different species), but all traps for the cameras
#plot(traps.tmp$posx,traps.tmp$posy)

#only expect to catch PEMA in the small trap grid, so need to specify that reduced trap file (P-Y)
traps.small <- traps.tmp[traps.tmp$Rlet > "H",] #removes traps from large grid
traps.large <- traps.tmp[traps.tmp$Rlet < "I",] #removes traps from small grid

traps.all <- traps.tmp #if using all traps
traps.cam <- traps.tmp[traps.tmp$cam==1,] #remove half of small grid for where no cameras

#check each species and capture locations
check <- merge(CR.data, traps.all, by.x="Loc", by.y="RC", all.x=T)

table(check[check$Species == "PEMA",][,c("ttype","Grid")])
table(check[check$Species == "TATO",][,c("ttype","Grid")])
table(check[check$Species == "GLSA",][,c("ttype","Grid")])

CR.data <- check[!(check$Species == "PEMA" & check$ttype == "tomahawk"),] #remove tomahawk traps for PEMA captures

#######
#formatting function
#######
format.fun<-function(type=c("cap"), data, traps, sp, t=NULL)
{ 
  #reformat capture data depending on data type
  #need an edf file
  if(type=="cap")
  {
    data$traps <- as.factor(paste(data$Stn,data$Stn.1, sep="")) #combine station letter and number
    #make edf, for indexing in package functions, need new trap numbers with no letters and in order match up correctly
    #length(levels(traps.tmp$RC))
    traps$trap <- 1:nrow(traps) #assign trap number to traps
    
    lookup <- data.frame(trap.LET=traps$RC, trap.NUM=traps$trap) #create match key
    
    data$trap <- with(lookup,trap.NUM[match(data$traps,trap.LET)]) #match up to assigned numbers
    data <- data[order(data$Grid),] #sort data by grid
    
    edf.data <- data.frame(session=data$Grid, individual=data$Combined_Tag, occasion=data$SHER_DAY, TrapID=data$trap) #encounter data frame
    tdf.data <- data.frame(TrapID=traps[order(traps$trap),]$trap, Xcoord=traps[order(traps$trap),]$posx,Ycoord=traps[order(traps$trap),]$posy) #trap data frame
    #,optional trap deployment matrix)) #replace with trap deployment file
    #if doing multi-strata, repeat the trap matrix above as the tdf list for length = number of strata
    table(edf.data[,c("session","occasion")])
    if(sp == "PEMA")
    {
      #use different values for session
      data$gnew <- data$Grid
      data[data$gnew == 7,]$gnew <- 6
      data[data$gnew == 8,]$gnew <- 7
      data[data$gnew == 9,]$gnew <- 8
      edf.data <- data.frame(session=data$gnew, individual=data$Combined_Tag, occasion=data$SHER_DAY, TrapID=data$trap) #encounter data frame
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
      scr.data <- data2oscr(edf=edf.data, 
                            sess.col=1, id.col=2, occ.col=3, trap.col=4, 
                            tdf=rep(list(tdf.data), 9),
                            K = c(8,8,8,4,8,6,8,7,4), #one per session
                            ntraps=rep(nrow(traps), 9), #one per session
                            remove.extracaps = FALSE) #counts for Poisson
    }
    sf <- scr.data$scrFrame
    sf
    # tiff(filename=paste(path.local, "Figures_SCR/", sp,"_SM-", grid, "_cr.tiff", sep=""), width=5, height=5, units="in", res=300, compression="lzw")
    # print(plot(sf))
    # dev.off()
    res <- 0.1*sf$mmdm
    ssDF <- make.ssDF(scrFrame=sf, buffer=round(3*sf$mmdm,0), res=res, cont.cov = F, fact.cov = F) #used res=20 since trap res was 10 at fine scale
    plot(ssDF)
    points(traps$posx,traps$posy)
    out.formatted <- list(sf=sf,ssDF=ssDF,res=res)
    #out.formatted <- list(sf=sf,ssDF=ssDF,res=res, S.Area=S.Area)
  }
  return(out.formatted)
}

######
#fit using oSCR (open spatial capture recapture)
######
oSCR.fun <- function(data)
{
  mod.oSCR <- oSCR.fit(list(D~session, p0~1, sig~1),data$sf,data$ssDF, encmod="P") #p0=encounter model, sig=spatial decay model, D=spatial point process model, density
  #these are for getting real of each grid separately
  #Dhat.oSCR <- exp(mod.oSCR$outStats$mle[3]) #this is in units of resolution (i.e., res)
  #Dhat.oSCR_m2 <- exp(mod.oSCR$outStats$mle[3])/data[[3]] #density/m^2
  #Dhat.oSCR_ha <- (exp(mod.oSCR$outStats$mle[3])/data[[3]])*10000 #density/ha
  #Nhat.oSCR <- exp(mod.oSCR$outStats$mle[3])*nrow(data$ssDF[[1]])
  out.oSCR <- list(mod=mod.oSCR, Dhat=Dhat.oSCR_ha,Nhat=Nhat.oSCR, nrows=nrow(data$ssDF[[1]]), area=data[[3]]) #you can add any other estiamtes and CIs here based on how you want to evaluate output
  #I have it returning density/ha so comparable across species as resolution is based on mmdm
  return(out.oSCR)
} 

###############
#fit using SCR_B
###############
SCR_B.fun <- function(data, n.chains, n.adapt, n.iter, n.burnin, n.cores, sp.cr)
  {
  #set up data
  y3d <- data$sf$caphist #extracting capture history array from oSCR scr frame
  S.Area <- nrow(data$ssDF[[1]])*data$res^2/10000 #this was in meters^2, now in ha
  mmax.g <- sapply(data$sf$caphist, function(x) nrow(x)*10)
  mmax <- max(mmax.g)
  #create 3D matrix: grid x ninds x traps
  yaug <- array(data=0, dim=c(length(y3d), mmax, ncol(data$sf$caphist[[1]]))) #dim(num sessions, mmax, num traps)
  for(g in 1:length(y3d))
  {
    yaug[g, 1:nrow(y3d[[g]]),] <- apply(y3d[[g]], 1:2, sum)
  }   
  
  locs <- cbind(data$ssDF[[1]]$X, data$ssDF[[1]]$Y)
  
  data.SCR_B <- list(n.sess=length(data$sf$caphist), mmax=mmax, J=ncol(data$sf$caphist[[1]]),
                     K=matrix(data=rep(data$sf$occasions, times=ncol(data$sf$caphist[[1]])), nrow=length(data$sf$caphist), ncol=ncol(data$sf$caphist[[1]])),
                     y=yaug, X=data$sf$traps[[1]], xlims=c(min(data$ssDF[[1]]$X), max(data$ssDF[[1]]$X)), 
                     ylims=c(min(data$ssDF[[1]]$Y), max(data$ssDF[[1]]$Y)), 
                     S.Area=S.Area)
  #initial values
  inits <- function(){list(psim=runif(length(data$sf$caphist)),
                           zm=array(data=1, dim=c(length(y3d), mmax)), #multi-session, start all m with zm=1
                           S=array(data=locs[sample(nrow(locs), size=mmax*length(y3d), replace=T),], dim=c(length(y3d), mmax, 2))
  )}
  
  #parameters to monitor
  params <- c("allgrids.mu.lam","sd.log.lam0","allgrids.mu.sigma","sd.log.sigma","lam0","sigma","psim", "N", "D", "S")
  
  save(data.SCR_B, inits, params, sp.cr, file=paste("RData/CR-SCR/SCR_randomeff_jagsargs_", sp.cr, ".RData", sep=""))
}

###########
#run this in parallel
###########
p <- function(sp.cr, traps) #g = grid, sp.cr = species name
{
  data.cr <- CR.data[which(CR.data$Species==sp.cr),]
  data.cr <- data.cr[data.cr$Loc %in% traps$RC,]
  if(nrow(data.cr) < 2)
  {
    write.table(paste("not enough", sp.cr, "individuals/captures to format grid", g, Sys.Date()), "test.txt", sep=",", append=T, col.names=F)
  }
  else
  {
    #remove data rows with no ear tags
    cr <-format.fun(type="cap", data=data.cr[!data.cr$Combined_Tag=="LR",], traps=traps, sp=sp.cr) #switch to traps = large if the tomahawks, switch traps= small if sherman
    
    out <- capture.output(cr)
    cat(out, file=paste("Output_SCR/", sp.cr, "_cr_", Sys.Date(), ".txt", sep=""), sep="\n", append=F)
    write.table(paste("Output_SCR/",sp.cr, "formatted grid for SCR ", Sys.Date()), "test.txt", sep=",", append=T, col.names=F)
    
    SCR_B <- SCR_B.fun(data=cr, n.chains=3, n.adapt = 500, n.iter=2500, n.burnin=500, n.cores=3, sp.cr=sp.cr)
  }
}

p(sp.cr="PEMA", traps=traps.small)

p(sp.cr="TATO", traps=traps.large)
p(sp.cr="GLSA", traps=traps.large)
