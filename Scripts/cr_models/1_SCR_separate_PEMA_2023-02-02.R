#updated 2019-09-09
#updated 2023-02-02

#SCR and SCC functions
#D. Morin 20190624

path.local <- "C:/Users/tosam/Documents/0_OSU/Dissertation/Manuscripts/SmallMammalCameras/SmammalCam_pub/"
# path.local <- "/raid1/home/fw/tosam/SmallMammalCams/"
# path.local <- "/nfs1/FW_HMSC/Levi_Lab/Tosa/Smammal_SCR/"
setwd(path.local)

######
#install.packages("devtools", repos='http://cran.us.r-project.org')
#install.packages("jagsUI", repos='http://cran.us.r-project.org')
#install.packages("MCMCglmm", repos='http://cran.us.r-project.org')

#library(devtools)
#install_github("jaroyle/oSCR")

######
library(oSCR)
library(jagsUI) #have to install jags directly on your machine
library(MCMCglmm)

require(plyr)

#20190620 Notes from MARIE on grids:

# Columns Rlet and Cnum in the SMcamera file represent trap stations.
# Columns Stn and Stn in the WOOPS live trap file represent trap stations.
# Rlet and Stn letters belong on the y axis. Letters A-H should correspond to positions: 280,240,200,160,120,80,40,0 OR  rev(seq(from=0, by=40, length.out = 8)); Letters P-Y should correspond to positions: 160,150,140,130,120,110,100,90,80,70 OR seq(from=160, by=-10, length.out=10)
# Cnum and stn numbers belong on the x axis. Numbers for stations A-H 1-8 correspond to positions: 0,40,80,120,160,200,240,280 OR seq(0, 280, by=40). Numbers for stations P-Y 1-10 correspond to positions: 160,170,180,190,200,210,220,230,240,250 OR seq(160, by=10, length.out=10)

#load and see what each of these files look like.
CR.data <- read.csv(paste(path.local,"Data_Raw/WOOPS_2017_proofed_2019-06-10.csv", sep=""), header=TRUE, sep=",")
CR.data$Loc <- paste(CR.data$Stn, CR.data$Stn.1, sep="")
#table(CR.data$Loc, CR.data$Grid)
#unique(CR.data$Loc)
CR.data <- CR.data[!CR.data$Loc=="NA",] #remove data rows with no trap location info
CR.data <- CR.data[!is.na(CR.data$SHER_DAY),] #remove data rows with no capture occasion number info

traps.tmp <- read.csv(paste(path.local, "Data_Raw/traplocs.txt", sep=""), sep=",") #this includes both a large grid and a small grid for the CR data 
#(different trap types catch different species), but all traps for the cameras
#plot(traps.tmp$posx,traps.tmp$posy)

#table(CR.data$Species)

#only expect to catch PEMA in the small trap grid, so need to specify that reduced trap file (P-Y)
traps.small <- traps.tmp[traps.tmp$Rlet > "H",] #removes traps from large grid
traps.large <- traps.tmp[traps.tmp$Rlet < "I",] #removes traps from small grid
traps.all <- traps.tmp #if using all traps
traps.cam <- traps.tmp[traps.tmp$cam==1,] #remove half of small grid for where no cameras

#formatting function
format.fun<-function(type=c("cap","cam"), data, traps, sp, grid, t=NULL)
{
  #reformat capture data depending on data type
  #need an edf file
  if(type=="cap"){
    data$traps <- as.factor(paste(data$Stn,data$Stn.1, sep="")) #combine station letter and number
    #make edf, for indexing in package functions, need new trap numbers with no letters and in order match up correctly
    #length(levels(traps.tmp$RC))
    #traps$trap <- as.integer(traps$RC) #reassign trap numbers in CR data based on traps.tmp$traps
    traps$trap <- 1:nrow(traps) #assign trap number to traps
    
    lookup <- data.frame(trap.LET=traps$RC, trap.NUM=traps$trap) #create match key
    
    data$trap <- with(lookup,trap.NUM[match(data$traps,trap.LET)]) #match up to assigned numbers
    edf.data<-data.frame(session=as.integer(1), individual=data$Combined_Tag, occasion=data$SHER_DAY, TrapID=data$trap) 
    #I only have data for one session, but if choose to process all sessions for multi-strata model change up session, but you'll want it to be a vector of integers 1:session to work best 
    tdf.data<-list(data.frame(TrapID=traps[order(traps$trap),]$trap, Xcoord=traps[order(traps$trap),]$posx,Ycoord=traps[order(traps$trap),]$posy))
    #,optional trap deployment matrix)) #replace with trap deployment file
    #if doing multi-strata, repeat the trap matrix above as the tdf list for length = number of strata
    scr.data <- data2oscr(edf=edf.data,sess.col=1,id.col=2,occ.col=3,trap.col=4,tdf=tdf.data,K = c(8),ntraps=c(nrow(traps)),
                          remove.extracaps = FALSE) #counts for Poisson
    sf <- scr.data$scrFrame
    sf
    tiff(filename=paste(path.local, "Figures/", sp,"_SM-", grid, "_cr.tiff", sep=""), width=5, height=5, units="in", res=300, compression="lzw")
    print(plot(sf))
    dev.off()
    ssDF <- make.ssDF(scrFrame=sf, buffer=round(2*sf$mmdm,0), res=0.1*sf$mmdm, cont.cov = F, fact.cov = F) #used res=20 since trap res was 10 at fine scale
    plot(ssDF)
    points(traps$posx,traps$posy)
    out.formatted<-list(sf=sf,ssDF=ssDF,res=0.1*sf$mmdm)
  }
  return(out.formatted)
}

######
#fit using oSCR (open spatial capture recapture)
oSCR.fun<-function(data){
  mod.oSCR<-oSCR.fit(list(D~1, p0~1, sig~1),data$sf,data$ssDF, encmod="P")
  Dhat.oSCR<-exp(mod.oSCR$outStats$mle[3]) #this is in units of resolution (i.e., res)
  Dhat.oSCR_m2<-exp(mod.oSCR$outStats$mle[3])/data[[3]] #density/m^2
  Dhat.oSCR_ha<-(exp(mod.oSCR$outStats$mle[3])/data[[3]])*10000 #density/ha
  Nhat.oSCR<-exp(mod.oSCR$outStats$mle[3])*nrow(data$ssDF[[1]])
  out.oSCR<-list(mod=mod.oSCR,Dhat=Dhat.oSCR_ha,Nhat=Nhat.oSCR)#you can add any other estiamtes and CIs here based on how you want to evaluate output
  #I have it returning density/ha so comparable across species as resolution is based on mmdm
  return(out.oSCR)
} 

###############
#fit using SCR_B
SCR_B.fun<-function(data, n.chains, n.adapt, n.iter, n.burnin, n.cores, g){
  sink("SCR.txt")
  cat("
        
        model{
        #priors
        lam0~dunif(0,5)
        sigma~dunif(0,100) #might need to adjust these based on what output looks like, I'm setting at approx 3*mmdm
        psim~dunif(0,1)
        
        #### model for marked individuals
        
        for (i in 1:mmax){ #
        zm[i] ~ dbern(psim) 
        S[i,1] ~ dunif(xlims[1], xlims[2])
        S[i,2] ~ dunif(ylims[1], ylims[2])
        
        for(j in 1:J){
        D2[i,j]<-(S[i,1]-X[j,1])^2 + (S[i,2]-X[j,2])^2
        lam[i,j]<-lam0*exp(-D2[i,j]/(2*sigma^2))
        y[i,j]~dpois(lam.effm[i,j]*K[j])  #model accumulated counts across K for marked individuals and add for whether a camera is functioning or not
        lam.effm[i,j]<-lam[i,j]*zm[i]
        }
        }
        
        N<-sum(zm[1:mmax])
        D<-N/S.Area
        
        } #end model description
        
        ",fill = TRUE)
  sink()
  
  #set up data
  y3d <- data$sf$caphist[[1]] #extracting capture history array from oSCR scr frame
  S.Area <- nrow(data$ssDF[[1]])*data$res^2/10000 #this was in meters, now in ha
  mmax <- nrow(data$sf$caphist[[1]])*10 #setting 5 times the number of individuals captured, 
  #but you can play with that if it looks like the posterior is bumping up against that max
  yaug <- matrix(0, nrow=mmax, ncol=ncol(data$sf$caphist[[1]])) #formatting 0s for data augmentation
  yaug[1:nrow(y3d), ] <- apply(y3d,1:2,sum )  #and adding back in real individuals
  
  data.SCR_B <- list(mmax=mmax, J=ncol(data$sf$caphist[[1]]), K=rep(data$sf$occasions, times=ncol(data$sf$caphist[[1]])),
                     y=yaug, X=data$sf$traps[[1]], xlims=c(min(data$ssDF[[1]]$X), max(data$ssDF[[1]]$X)), 
                     ylims=c(min(data$ssDF[[1]]$Y), max(data$ssDF[[1]]$Y)), 
                     S.Area=S.Area)
  #initial values
  inits <- function(){list(lam0=runif(1), psim=runif(1),
                           zm=rep(1,data.SCR_B$mmax), 
                           S=cbind(data$ssDF[[1]]$X,data$ssDF[[1]]$Y)[sample(nrow(cbind(data$ssDF[[1]]$X,data$ssDF[[1]]$Y)),
                                                                             size=data.SCR_B$mmax,replace=TRUE),])} #sampling with replacement
  #parameters to monitor
  params <- c("lam0","sigma","psim", "N", "D")
  
  save.image(paste(path.local, "RData/SCR_separate_", sp.cr, "_", g, ".RData", sep=""))
  
  #test<-jags.model(file="SCR.txt",data=data.SCR_B, inits=inits, n.chain=n.chain, n.adapt=n.adapt)
  
  mod.SCR_B<-jags(data=data.SCR_B, inits=inits, parameters.to.save = params, model.file="SCR.txt", n.chains=n.chains,
                  n.adapt=n.adapt, n.iter=n.iter, n.burnin=n.burnin, n.thin=1,
                  parallel = TRUE, n.cores = n.cores)
  D.mode<-posterior.mode(as.mcmc(mod.SCR_B[[1]]$D)) #*10000 #/ha
  N.mode<-posterior.mode(as.mcmc(mod.SCR_B[[1]]$N))
  
  out.SCR_B<-list(mod=mod.SCR_B, Dhat=D.mode, Nhat=N.mode, SArea=S.Area)
  return(out.SCR_B)
}

##########
table(CR.data$Species, CR.data$Grid)
########

#run this in parallel
#for(g in unique(CR.data$Grid))
p <- function(g)
{
  print(g)
  data.cr <- CR.data[which(CR.data$Species==sp.cr & CR.data$Grid==g),]
  data.cr <- data.cr[data.cr$Loc %in% traps$RC,]
  
  if(nrow(data.cr) < 2)
  {
    write.table(paste("not enough", sp.cr, "individuals/captures to format grid", g, Sys.Date()), "test.txt", sep=",", append=T, col.names=F)
  }
  else
  {
    #####
    #remove data rows with no ear tags
    cr <- format.fun(type="cap", data=data.cr[!data.cr$Combined_Tag=="LR",], traps=traps, sp=sp.cr, grid=paste(0, g, sep="")) #switch to traps = large if the tomahawks, switch traps= small if sherman
    out <- capture.output(cr)
    cat(out, file=paste("Output_SCR/", sp.cr, "_SM-0", g, "_cr_", Sys.Date(), ".txt", sep=""), sep="\n", append=F)
    write.table(paste(sp.cr, "formatted grid for SCR ", g, Sys.Date()), "test.txt", sep=",", append=T, col.names=F)
    
    #Spatial Capture Recapture Max Likelihood Framework
    #oSCR<-oSCR.fun(data=cr)
    #out <- capture.output(oSCR)
    #cat(out, file=paste("Output/", sp.cr, "_SM-0", g, "_oSCR_", Sys.Date(), ".txt", sep=""), sep="\n", append=F)
    
    #Spatial Capture Recapture, Bayesian Framework
    SCR_B <- SCR_B.fun(data=cr, n.chains=3, n.adapt = 500, n.iter=2500, n.burnin=500, n.cores=3, g=g)  
    out <- capture.output(SCR_B)
    cat(out, file=paste("Output/", sp.cr, "_SM-0", g, "_SCR_B_", Sys.Date(), ".txt", sep=""), sep="\n", append=F)
    
    tiff(filename=paste("Figures/Traceplots_", sp.cr, "_SM-0", g, "_SCR-B.cr_", Sys.Date(), ".tiff", sep=""), width=5, height=8, units="in", res=300, compression="lzw")
    print(traceplot(SCR_B$mod$samples))
    dev.off()
  }
}

require(parallel)
no_cores <- detectCores() - 1
# no_cores <- 9
cl <- makeCluster(no_cores)
clusterEvalQ(cl, library(oSCR))
clusterEvalQ(cl, library(jagsUI))
clusterEvalQ(cl, library(MCMCglmm))
clusterEvalQ(cl, library(plyr))

sp.cr <- "PEMA"
sp.cam <- "DEER MOUSE"
traps=traps.small
#mclapply(1:5, p, mc.cores=5)

clusterExport(cl, c("CR.data","sp.cr","sp.cam","format.fun","traps","traps.cam","path.local","p","oSCR.fun","SCR_B.fun"))
glsa_scr <- parLapply(cl=cl, c(1:9), fun=p)

stopCluster(cl)
