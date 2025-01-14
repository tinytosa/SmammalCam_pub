#Nmixture models with jags

#######################
#load packages
require(dplyr)
require(plyr)
require(reshape2)
require(ggplot2)

#######################
#load data
path.local <- getwd()
setwd(path.local)

#import small mammal data
smdata <- read.table("Data_Raw/SMcamera_tags_daynum_2020-05-04.txt", sep=",") #only includes data for deer mice, flying squirrels, townsends chipmunk, and voles
smdata$DT <- as.POSIXct(smdata$DT, format="%Y-%m-%d %H:%M:%S", tz="GMT")

traps.tmp <- read.csv("Data_Raw/traplocs.txt", sep=",") #this includes both a large grid and a small grid for the CR data 

#as.integer for characters stopped working
traps.small <- traps.tmp[traps.tmp$Rlet > "H",]
traps.large <- traps.tmp[traps.tmp$Rlet < "I",]

traps.all <- traps.tmp #if using all traps
traps.cam <- traps.tmp[traps.tmp$cam == 1,] #remove half of small grid for where no cameras

traps.both <- merge(traps.large, traps.small, by=c("posx","posy"), suffixes=c(".t",".s")) #find RC of tomahawk traps that correspond to sherman grid
traps.cam.small <- traps.cam[traps.cam$RC %in% traps.both$RC.t | traps.cam$RC %in% traps.small$RC,] #nested cameras along sherman grid
traps.cam.large <- traps.cam[traps.cam$ttype == "tomahawk",]

traps.t <- data.frame(traps.cam.large %>% slice(rep(1:n(), times = 8)), grid=rep(1:8, each=64))
traps.t$Directory <- paste(traps.t$grid, traps.t$RC, sep="")

traps.s <- data.frame(traps.cam.small %>% slice(rep(1:n(), times = 8)), grid=rep(1:8, each=25))
traps.s$Directory <- paste(traps.s$grid, traps.s$RC, sep="")

########
#only use shermans for PEMA
#only use tomahawks for TATO and GLSA


# require(xlsx)
#import cmr abundance estimates
#huggins estimates
#cmr <- read.xlsx("SmallMammal_Huggins_Estimates_fromWeldy_updated.xlsx", sheetName="Abundances")
#cmr$sp <- tolower(cmr$Sp)
#cmr <- cmr[cmr$ModelName == "Huggins Null",]

#SCR_B estimates
# cmr <- read.table("Output_SCR_B_2020-07-03.txt", sep=",", header=T)
# cmr$sp <- tolower(cmr$sp)
# cmr$Grid <- substr(cmr$grid, 5, 5)

#merge with oSCR results
# cmr <- read.table("Output_oSCR_allgrids_2020-08-26.txt", sep=",", header=T)
# cmr$g <- cmr$grid
# cmr$grid <- substr(cmr$grid, 5, 6)
# cmr$sp <- factor(cmr$sp, levels=c("PEMA","TATO","GLSA"))
# cmr$scaledscrn <- cmr$scaledn

# cmr <- read.table("Output_SCR_randomeff_2020-08-12.txt", sep=",", header=T)
# cmr$sp <- tolower(cmr$sp)

keys <- data.frame(key = unique(smdata$Keywords))
keys$abbr <- c("tato","pema","glsa","myo")
keys$a <- c("t","p","g","m")
#keys$lambdamax <- c(100, 300, 50, 50) #for initial lambda runif max parameter

# remove sherman grid cameras?
# smdata.t <- smdata[!smdata$Rlet %in% c("P","R","T","V","X"),] #only tomahawk grid cameras

########
#format data so each row is a site, each column is a rep, separate dataframe for each grid

#key <- "DEER MOUSE"
#grid <- 1

for(t in c(0, 15, 60, 1440))
{
  print(t)
  for(key in keys$key)
  {
    print(key)
    k <- keys[keys$key == key,]$abbr #get abbreviation
    if(k == "pema")
    {
      traps <- traps.s
    }
    else
    {
      traps <- traps.t
    }
    smdata.k <- smdata[smdata$Directory %in% traps$Directory,]
    assign(k, smdata.k[smdata.k$Keywords == key,]) #subset data for keyword, pema
    for(grid in 1:8)
    {
      print(grid)
      #assign(paste(k,grid, sep="_"), eval(parse(text=k))[eval(parse(text=k))$grid==grid,]) #subset for grid, pema_1
      g <- eval(parse(text=k))[eval(parse(text=k))$grid==grid,] #subset for grid
      
      #sort by date, time, and location
      g <- g[order(g$DT),]
      g <- g[order(g$Directory),]
      
      #consolidate photos within t minutes
      g$timediff <- c(0,round((g$DT[-1] - g$DT[-nrow(g)])/60, digits=0)) #need to do this for each species
      
      n <- 1
      g$consolidate <- sapply(g$timediff, function(x) {if(x > t | x < 0)
      {n <<- n+1} 
        n})
      
      #consolidate if keywords and consolidate numbers are the same
      g <- ddply(g, .(Directory, grid, Loc, Rlet, Cnum, Keywords, num, consolidate, nocturnal), summarise, numphoto=length(DT), startDT=min(DT))
      
      #create dataframe with each loc as row and each rep as column
      #assign(paste(k,grid,"l", sep="_"), ddply(eval(parse(text=paste(k,grid, sep="_"))), .(Directory, num), nrow)) #long format
      l <- ddply(g, .(Directory, num), nrow) #long format
      l <- merge(l, traps[traps$grid==grid,], by="Directory", all.y=T) #make sure data frame includes all functional cameras
      l[is.na(l)] <- 0
      #l <- l[!l$num == 0,]
      #assign(paste(k,grid,"l", sep="_"), l)
      #names: Directory, num, V1
      w <- dcast(l[,!names(l) %in% "grid"], Directory ~ num, value.var="V1", fill=0)
      w <- w[,!names(w) %in% 0]
      assign(paste(k,grid,"w_t",t, sep="_"), w) #wide format
      # write.table(file=paste("Data_Raw/data_perdaynum/",k,"_",grid,"_t",t,".txt", sep=""), w, sep=",", row.names=F) #write out as file, only need to do this once
    }
  } 
}

#######
#0. load data
# pema_1_w <- read.table("Data_Raw/data_perdaynum/pema_1_t0.txt", sep=",", header=T)
# glsa_1_w <- read.table("Data_Raw/data_perdaynum/glsa_1_t0.txt", sep=",", header=T)

#list files
files <- dir("Data_Raw/data_perdaynum", pattern=".txt")

#load all data
for(f in files)
{
  assign(gsub(f, pattern=".txt", replacement="_w"), read.table(paste("Data_Raw/data_perdaynum/",f, sep=""), sep=",", header=T))
}

#######
#1. basic model
#######

sink("Models/base.txt")
cat("
model{
  for (i in 1:Nsites) #for each site
  {
    N[i]~dpois(lambda) #state model

    #likelihood
    for (j in 1:Nreps)
    {
      y[i,j]~dbinom(p,N[i]) #observation model
    }
  }
  
  #priors
  #p~dbeta(1,1)
  p~dunif(0,1) #uniform
  lambda~dgamma(0.001,0.001) #underlying density

  sumN <- sum(N[])/Nsites #calculate average
}
", fill=TRUE)
sink()

############
#1a. apply basic model to glsa, tato, pema, myo data
############

require(jagsUI)
require(lattice)
require(coda)

dname <- ls(pattern="_w")
#d <- dname[10] #myo_0
#d <- dname[2] #glsa 1
#dname[19:27] #pema

for(d in dname)
{
  print(d)
  data_m <- as.matrix(eval(parse(text=d))[,!names(eval(parse(text=d))) %in% "Directory"]) #matrix format
  #data_m <- data_m[,-1] #remove day 1?
  Nmix1.data <- list(y=data_m, Nsites=nrow(data_m), Nreps=ncol(data_m))
  Ninit <- apply(data_m, 1, max) #get max count at each site
  
  Init.fun <- function() #calculate initial values for parameters
  {
    list(N=Ninit,"lambda"=runif(1,1,500),"p"=runif(1,0.1,0.9)) #if N not given, one may be chosen where obs>N*(p=1) and 5~dbino(3,1) is error
  }
  #params=c("lambda","p")
  params=c("lambda","p","N", "sumN")
  
  jags.out <- jags(model.file="Models/base.txt", data=Nmix1.data, inits=Init.fun,
                   parameters.to.save=params, n.chains=3, n.iter=20000, n.burnin=10000, n.thin=5, n.adapt=20000)
  print(jags.out)
  
  # write.table(jags.out$summary, file=paste("Output_Nmix/outjagsUI_base_params_", gsub(d, pattern="_w", replacement=""), "_", Sys.Date(), ".txt", sep=""), sep=",")
  
  #traceplot(jags.out) #plots traceplots one after another, need to hit enter
  #densityplot(jags.out)
  
  jagsout.mcmc <- as.mcmc.list(jags.out$samples)
  p <- xyplot(jagsout.mcmc[,c("lambda", "p","sumN")]) #plots traceplots in one figure
  # dp <- densityplot(jagsout.mcmc[,c("lambda", "p","sumN")]) #density plot of posterior
  
  # tiff(filename=paste("Figures_Nmix/outjags_trace_",gsub(d, pattern="_w", replacement=""), "_", Sys.Date(), ".tiff", sep=""), height=6, width=6, res=300, units="in", compression="lzw")
  # print(p)
  # dev.off()
  
  # tiff(filename=paste("Figures_Nmix/outjags_densityplot_",gsub(d, pattern="_w", replacement=""), "_", Sys.Date(), ".tiff", sep=""), height=6, width=6, res=300, units="in", compression="lzw")
  # print(dp)
  # dev.off()
}

############
#1b. look at nmixture model output
#convergence??
############
require(tidyr)

setwd("C:/Users/tosam/Documents/0_OSU/Dissertation/Data/SmallMammalCams/")

o.nmix <- data.frame(fullfilename=dir("Data_Model_Output/Output_Nmix/alldays_base", pattern=".txt", full.names = T))
o.nmix <- separate(o.nmix, col="fullfilename", into=c("Data_Model_Output","outfolder","model","filename"), sep="/", remove=F)
o.nmix <- separate(o.nmix, col="filename", into=c("outjags","det","params","sp","grid","t","date"), sep="_", remove=F)
o.nmix$date <- gsub(o.nmix$date, pattern=".txt", replacement="")

out.all <- data.frame()

for(o in o.nmix$fullfilename)
{
  print(o)
  
  out <- read.table(as.character(o), sep=",")
  out$date <- o.nmix[o.nmix$fullfilename==o,]$date
  out$sp <- o.nmix[o.nmix$fullfilename==o,]$sp
  out$grid <- o.nmix[o.nmix$fullfilename==o,]$grid
  out$param <- row.names(out)
  
  out.all <- rbind(out.all, out)
}

out.all$converged <- "not converged"
out.all[out.all$Rhat < 1.1,]$converged <- "converged"

unique(out.all[out.all$Rhat > 1.1,][,c("grid","sp")])

#########
#1c.plot nmixture parameters against CMR estimates
#########
lambda <- out.all[out.all$param == "lambda",]

# tiff(filename="Figures_Nmix/nmix_summary_scr_alldays.tiff", width=6, height=8, units="in", res=300, compression="lzw")
ggplot(data=lambda[!lambda$grid == 0,]) + 
  geom_point(aes(x=mean, y=grid, group=sp, col=converged), size=3) + 
  geom_segment(aes(x=X2.5., xend=X97.5., y=grid, yend=grid, col=converged), lwd=1) +
  # geom_point(data=cmr[!cmr$Grid ==9,], aes(x=nmode, y=as.numeric(Grid)+0.25, group=sp), size=3, color="steelblue", pch=17) + #for SCR_B data
  # geom_errorbar(data=cmr[!cmr$Grid ==9,], aes(xmin=n95low, xmax=n95high, y=as.numeric(Grid)+0.25, group=sp), color="steelblue", width=0.1) +
  # geom_point(data=cmr[!cmr$grid ==9,], aes(x=n, y=as.numeric(grid)+0.25, group=sp), size=3, color="steelblue", pch=17) + #for oSCR_all data
  # geom_errorbar(data=cmr[!cmr$grid ==9,], aes(xmin=n_lwr, xmax=n_upr, y=as.numeric(grid)+0.25, group=sp), color="steelblue", width=0.1) +
  
  scale_x_log10() + 
  scale_color_manual(values=c("black","grey80"), name="") +
  facet_wrap(~sp, scales="free", ncol=1) + theme_bw() + theme(legend.position="top")
# dev.off()

#######
#1extra. plot abundances
#######
require(tidyr)
require(ggplot2)
require(ggpubr)

#plot detected values
#g <- 1
#sp <- "glsa"
for(g in 1:8)
{
  for(sp in c("pema","glsa","tato"))
  {
    #plot detection values
    og_data <- data.frame(Directory=eval(parse(text=paste(sp, g, "t0_w", sep="_")))[,1], total=rowSums(eval(parse(text=paste(sp, g, "t0_w", sep="_")))[,-1]))
    og_data$gridloc <- substr(og_data$Directory, 2, 3)
    og_data <- merge(og_data, traps.tmp, by.x="gridloc", by.y="RC", all.x=T)
    
    p.det <- ggplot(data=og_data) + geom_tile(aes(posx, y=posy, fill=total)) + scale_fill_distiller(palette="Spectral") + theme_bw(base_size = 20) + coord_equal() + theme(legend.position="top")
    # ggsave(p.det, filename=paste("Figures/SM_cams_detections_", sp, "_g", g, ".tiff", sep=""), height=6, width=6, units="in", dpi=300, compression="lzw")
  }
}

#plot predicted n from model
for(g in 1:8)
{
  for(sp in c("pema","glsa","tato"))
  {
    #plot predicted n values
    N.vals <- out.all[grep(out.all$param, pattern="N"),]
    N.g <- N.vals[N.vals$grid == g & N.vals$sp == sp,]
    N.g <- data.frame(loc=eval(parse(text=paste(sp, g, "t0_w", sep="_")))[,"Directory"], n=N.g$mean)
    N.g$gridloc <- substr(N.g$loc, 2, 3)
    N.g <- merge(N.g, traplocs, by.x="gridloc", by.y="RC", all.x=T)
    p.N <- ggplot(data=N.g) + geom_tile(aes(x=posx, y=posy, fill=n)) + 
      scale_fill_distiller(palette="Spectral") +
      theme_bw(base_size=10) + coord_equal() + theme(legend.position="top") + ggtitle(toupper(sp))
    ggsave(p.N, filename=paste("Figures_Nmix/jagsUI_n_", sp, "_g",g, ".tiff", sep=""), height=6, width=6, units="in", dpi=300, compression="lzw")
    assign(sp, p.N)
  }
  all.N <- ggarrange(pema, tato, glsa, nrow=1)
  ggsave(all.N, filename=paste("Figures_Nmix/jagsUI_n_all_g",g, ".tiff", sep=""), height=6, width=10, units="in", dpi=300, compression="lzw")
}

ggarrange(p.det, p.N)

#########
#1d. plot ratio of nmixture to CMR estimates
#########
t=0
lambda$sp <- toupper(lambda$sp)
#ratio <- merge(lambda, cmr, by.x=c("grid","sp"), by.y=c("Grid","sp"), all.x=T)
ratio <- merge(lambda, cmr, by.x=c("grid","sp"), by.y=c("grid","sp"), all.x=T)
ratio <- ratio[ratio$sp != "myo" & ratio$grid != 0,] #remove rows where no SCR data

#scale values for each species
ratio$scaled <- 0
#ratio$scaledN <- 0 #already scaled scr_n
for(sp in unique(ratio$sp))
{
  z <- ratio[ratio$sp == sp,]
  if(nrow(z) == 0)
  {next}
  ratio[ratio$sp == sp,]$scaled <- scale(z$mean)
  #ratio[ratio$sp == sp,]$scaledN <- scale(z$N)
  
  lm.scaled <- lm(data=ratio[ratio$sp == sp,], scaled ~ scaledscrn)
  sink(paste("Output_Nmix/nmix_scaled_lm_", sp, "_base_t", t, "_SCR.txt", sep=""))
  print(summary(lm.scaled))
  sink()
}

ratio$sp <- factor(ratio$sp,levels=c("pema","tato","glsa"))

#plot raw numbers
# r0 <- ggplot(data=ratio) + 
#   geom_point(aes(x=nmode, y=mean, group=sp, col=converged), size=3) +
#   stat_smooth(aes(x=nmode, y=mean), method="lm", formula= y~x, col="black") + 
#   scale_color_manual(values=c("black","grey80"), name="") +
#   facet_wrap(~sp, scale="free", nrow=1) +
#   theme_bw() + theme(legend.position="top")
# r0
# ggsave(r0, filename="Figures_Nmix/nmix_comparison_all.tiff", height=6, width=8, units="in", dpi=300, compression="lzw")

r00 <- ggplot(data=ratio) +
  geom_point(aes(x=nmode, y=scaled, group=sp), pch=1, col="grey80") +
  geom_point(aes(x=nmode, y=scaled, group=sp, col=converged), size=3) +
  stat_smooth(aes(x=nmode, y=scaledscrn), method="lm", formula= y~x, col="grey80") + 
  stat_smooth(aes(x=nmode, y=scaled), method="lm", formula= y~x, col="black") + 
  scale_color_manual(values=c("black","grey80"), name="") +
  facet_wrap(~sp, scale="free", nrow=1) +
  theme_bw() + theme(legend.position="top")
r00
ggsave(r00, filename="Figures_Nmix/nmix_scr_comparisonscaled_base_all.tiff", height=6, width=8, units="in", dpi=300, compression="lzw")

##############
############
#2. change detection probability by site
############

#logit.p[i] <- a0
#p[i]<- exp(logit.p[i])/(exp(logit.p[i])+1) #inverse logit


sink("Models/psite.txt")
cat("
model
{
  for (i in 1:Nsites) #for each site
  {
    p[i]~dunif(0,1)

    N[i]~dpois(lambda) #state model
    
    #likelihood
    for (j in 1:Nreps)
    {
      y[i,j]~dbinom(p[i], N[i]) #observation model
    }
  }
  
  #priors
  lambda~dgamma(0.001,0.001)
  
  sumN <- sum(N[])/Nsites #calculate average
}
", fill=TRUE)
sink()

########
#2b. run model2, change detection probability by site
########
dname <- ls(pattern="_w")
dname <- dname[!dname %in% grep(dname, pattern="myo", value=T)]
#dname <- dname[c(2,11,20,29)]
#d <- dname[10] #myo_0
#d <- dname[2] #glsa 1

for(d in dname)
{
  print(d)
  data_m <- as.matrix(eval(parse(text=d))[,!names(eval(parse(text=d))) %in% "Directory"]) #matrix format
  Nmix2.data <- list(y=data_m, Nsites=nrow(data_m), Nreps=ncol(data_m))
  
  Ninit2 <- apply(data_m,1,max) #get max count at each site
  
  Init2.fun<-function()
  {
    list(N=Ninit2, "lambda"=runif(1,1,500))
  }
  params=c("lambda", "p", "N", "sumN")
  
  jags.out <- jags(model.file="psite.txt", data=Nmix2.data, inits=Init2.fun,
                   parameters.to.save=params, n.chains=3, n.iter=15000, n.burnin=5000, n.thin=5, n.adapt=20000)
  print(jags.out)
  write.table(jags.out$summary, file=paste("Output_Nmix/outjagsUI_psite_params_", gsub(d, pattern="_w", replacement=""), "_", Sys.Date(), ".txt", sep=""), sep=",")
  
  #traceplot(jags.out$samples) #plots traceplots one after another, need to hit enter
  
  jagsout.mcmc <- as.mcmc.list(jags.out$samples)
  p <- xyplot(jagsout.mcmc[,c("lambda","p[1]","N[1]","sumN")]) #plots traceplots in one figure
  # dp <- densityplot(jagsout.mcmc[,c("lambda","p[1]","N[1]")]) #density plot of posterior
  
  tiff(filename=paste("Figures_Nmix/outjags_psite_trace_",gsub(d, pattern="_w", replacement=""), "_", Sys.Date(), ".tiff", sep=""), height=6, width=6, res=300, units="in", compression="lzw")
  print(p)
  dev.off()
  
  # tiff(filename=paste("Figures_Nmix/outjags_psite_densityplot_",gsub(d, pattern="_w", replacement=""), "_", Sys.Date(), ".tiff", sep=""), height=6, width=6, res=300, units="in", compression="lzw")
  # print(dp)
  # dev.off()
}

############
############
#3. detection is not equal across days
######
p1 <- ggplot(pema_1_l) + geom_path(aes(x=num, y=V1, group=Directory), lwd=1) + facet_wrap(~Rlet, ncol=2) + theme_bw()
ggsave(p1, filename="detectionsperday_pema_1.tiff", height=8, width=4, units="in", dpi=300, compression="lzw")

p2 <- ggplot(glsa_1_l) + geom_path(aes(x=num, y=V1, group=Directory), lwd=1) + facet_wrap(~Rlet, ncol=2) + theme_bw()
ggsave(p2, filename="detectionsperday_glsa_1.tiff", height=8, width=4, units="in", dpi=300, compression="lzw")

p3 <- ggplot(tato_1_l) + geom_path(aes(x=num, y=V1, group=Directory), lwd=1) + facet_wrap(~Rlet, ncol=2) + theme_bw()
ggsave(p3, filename="detectionsperday_tato_1.tiff", height=8, width=4, units="in", dpi=300, compression="lzw")

p4 <- ggplot(myo_1_l) + geom_path(aes(x=num, y=V1, group=Directory), lwd=1) + facet_wrap(~Rlet, ncol=2) + theme_bw()
ggsave(p4, filename="detectionsperday_myo_1.tiff", height=8, width=4, units="in", dpi=300, compression="lzw")

############
#3a. detection by site and decay over time
########
sink("Models/pdecay.txt")
cat("
model
{

  for (i in 1:Nsites) #for each site
  {
    det[i]~dunif(0,1)
    N[i]~dpois(lambda) #state model
    
    #likelihood
    for (j in 1:Nreps) #for each obs
    {
      p[i,j] <- det[i]*exp(d0*(j-1))
      y[i,j]~dbinom(p[i,j], N[i])
    }
  }
  
  #priors
  d0~dunif(-10,0)
  lambda~dgamma(0.001,0.001)
  
  sumN <- sum(N[])/Nsites #calculate average, derived parameters
}
", fill=TRUE)
sink()


#########
#3b. run jagsUI models for pdecay
#########
#d <- "tato_1_w"
dname <- ls(pattern="_w")
#dname <- dname[c(2,11,20,29)]
#d <- dname[10] #myo_0
#d <- dname[2] #glsa 1
#dname[19:27] #pema

for(d in dname)
{
  print(d)
  data_m <- as.matrix(eval(parse(text=d))[,!names(eval(parse(text=d))) %in% "Directory"]) #matrix format
  Nmix3.data <- list(y=data_m, Nsites=nrow(data_m), Nreps=ncol(data_m))
  
  Ninit3 <- apply(data_m,1,max)#get max count at each site
  
  Init3.fun<-function()
  {
    #list(N=Ninit3, "a0"=rnorm(1,0,1), "b0"=rnorm(1,0,1), "d0"=rnorm(1,0,1), "lambda"=runif(1,1,500))#if N not given, one may be chosen where obs>N*(p=1) and 5~dbino(3,1) is error
    #list(N=Ninit3, "a0"=1, "b0"=1, "d0"=-1, "lambda"=runif(1,1,500)) #if N not given, one may be chosen where obs>N*(p=1) and 5~dbino(3,1) is error
    list(N=Ninit3, "d0"=-1, "lambda"=runif(1,1,500))
  }
  #params=c("lambda","a0","b0","d0")
  params=c("lambda","d0","det","N","sumN")
  
  jags.out <- jags(model.file="pdecay.txt", data=Nmix3.data, inits=Init3.fun, 
                   parameters.to.save=params, n.chains=3, n.iter=15000, n.burnin=5000, n.thin=5, n.adapt=20000)
  print(jags.out)
  write.table(jags.out$summary, file=paste("Output_Nmix/outjagsUI_pdecay_params_", gsub(d, pattern="_w", replacement=""), "_", Sys.Date(), ".txt", sep=""), sep=",")
  
  #traceplot(jags.out$samples) #plots traceplots one after another, need to hit enter
  
  jagsout.mcmc <- as.mcmc.list(jags.out$samples)
  p <- xyplot(jagsout.mcmc[,c("lambda","d0","sumN")]) #plots traceplots in one figure
  # dp <- densityplot(jagsout.mcmc[,c("lambda","d0","sumN")]) #density plot of posterior
  
  tiff(filename=paste("Figures_Nmix/outjags_pdecay_trace_",gsub(d, pattern="_w", replacement=""), "_", Sys.Date(), ".tiff", sep=""), height=6, width=6, res=300, units="in", compression="lzw")
  print(p)
  dev.off()
  
  # tiff(filename=paste("Figures_Nmix/outjags_pdecay_densityplot_",gsub(d, pattern="_w", replacement=""), "_", Sys.Date(), ".tiff", sep=""), height=6, width=6, res=300, units="in", compression="lzw")
  # print(dp)
  # dev.off()
  
  #beta distribution shape
  #pbeta <- rbeta(100, jags.out$BUGSoutput$mean$a0, jags.out$BUGSoutput$mean$b0)
  #ggplot() + geom_histogram(aes(x=pbeta), binwidth=0.025) + theme_bw()
}

############
#3c. look at nmixture model outputs of pdecay
#convergence??
############
require(tidyr)
# o.nmix <- data.frame(fullfilename=dir("Output_Nmix/UI_pdecay_alldays_2020-06-26/", pattern=".txt", full.names = T))
o.nmix <- data.frame(fullfilename=dir("Data_Model_Output/Output_Nmix/alldays_pdecay/", pattern=".txt", full.names = T))
o.nmix <- separate(o.nmix, col="fullfilename", into=c("outfolder","output_nmix","model","filename"), sep="/", remove=F)
o.nmix <- separate(o.nmix, col="filename", into=c("outjags","pdecay","params","sp","grid","t","date"), sep="_", remove=F)
o.nmix$date <- gsub(o.nmix$date, pattern=".txt", replacement="")

out.all <- data.frame()

for(o in o.nmix$fullfilename)
{
  print(o)
  
  out <- read.table(o, sep=",")
  out$date <- o.nmix[o.nmix$fullfilename==o,]$date
  out$sp <- o.nmix[o.nmix$fullfilename==o,]$sp
  out$grid <- o.nmix[o.nmix$fullfilename==o,]$grid
  out$param <- row.names(out)
  
  out.all <- rbind(out.all, out)
}

out.all$converged <- "not converged"
out.all[out.all$Rhat < 1.1,]$converged <- "converged"

unique(out.all[out.all$Rhat > 1.1,][,c("grid","sp")])

###########
#plot abundances
###########
require(tidyr)
require(ggplot2)
require(ggpubr)

#plot predicted n from model

for(g in 0:8)
{
  for(sp in c("pema","glsa","tato"))
  {
    #plot predicted n values
    N.vals <- out.all[grep(out.all$param, pattern="N"),]
    N.g <- N.vals[N.vals$grid == g & N.vals$sp == sp,]
    N.g <- data.frame(loc=eval(parse(text=paste(sp, g, "w", sep="_")))[,"Directory"], n=N.g$mean)
    N.g$gridloc <- substr(N.g$loc, 2, 3)
    N.g <- merge(N.g, traplocs, by.x="gridloc", by.y="RC", all.x=T)
    p.N <- ggplot(data=N.g) + geom_tile(aes(x=posx, y=posy, fill=n)) + 
      scale_fill_distiller(palette="Spectral") +
      theme_bw(base_size=10) + coord_equal() + theme(legend.position="top") + ggtitle(toupper(sp))
    ggsave(p.N, filename=paste("Figures_Nmix/jagsUI_n_pdecay", sp, "_g",g, ".tiff", sep=""), height=6, width=6, units="in", dpi=300, compression="lzw")
    assign(sp, p.N)
  }
  all.N <- ggarrange(pema, tato, glsa, nrow=1)
  ggsave(all.N, filename=paste("Figures_Nmix/jagsUI_n_pdecay_all_g",g, ".tiff", sep=""), height=6, width=10, units="in", dpi=300, compression="lzw")
}

##########
det.all <- out.all[grep(out.all$param, pattern="det"),]

ggplot(det.all, aes(x=mean, y=param, group=sp, col=sp)) + geom_point() + facet_wrap(~grid) + theme_bw()

##########
#3d. plot nmixture parameters against CMR estimates
##########
require(ggplot2)
require(xlsx)
# cmr <- read.xlsx("SmallMammal_Huggins_Estimates_fromWeldy_updated.xlsx", sheetName="Abundances")
# cmr$sp <- tolower(cmr$Sp)
# cmr <- cmr[cmr$ModelName == "Huggins Null",]

lambda <- out.all[out.all$param == "lambda",]
#merge(lambda, cmr, )

tiff(filename="Figures_Nmix/nmix_summary_pdecay_alldays.tiff", width=6, height=8, units="in", res=300, compression="lzw")
ggplot(data=lambda[!lambda$grid == 0,]) + 
  geom_point(aes(x=mean, y=grid, group=sp, col=converged), size=3) + 
  geom_segment(aes(x=X2.5., xend=X97.5., y=grid, yend=grid, col=converged), lwd=1) +
  geom_point(data=cmr[cmr$ModelName=="Huggins Null" & !cmr$Grid == 9,], aes(x=N, y=Grid+0.25, group=sp), size=3, color="steelblue", pch=17) +
  scale_x_log10() + 
  scale_color_manual(values=c("black","grey80"), name="") +
  facet_wrap(~sp, scales="free", ncol=1) + theme_bw() + theme(legend.position="top")
dev.off()

#########
#3e. plot ratio of nmixture to CMR estimates
#########
ratio3 <- merge(lambda, cmr, by.x=c("grid","sp"), by.y=c("Grid","sp"), all.x=T)
ratio3 <- ratio3[!ratio3$sp %in% c("myo") & !ratio3$grid %in% c(0),]

#scale values for each species
ratio3$scaled <- 0
ratio3$scaledN <- 0

for(y in unique(ratio3$sp))
{
  z <- ratio3[ratio3$sp == y,]
  if(nrow(z) == 0)
  {next}
  ratio3[ratio3$sp == y,]$scaled <- scale(z$mean)
  ratio3[ratio3$sp == y,]$scaledN <- scale(z$N)
}

ratio3$sp <- factor(ratio3$sp,levels=c("pema","tato","glsa"))

r3 <- ggplot(data=ratio3) + 
  geom_point(aes(x=N, y=mean, group=sp, col=converged), size=3) +
  stat_smooth(aes(x=N, y=mean), method="lm", formula= y~x) + 
  scale_color_manual(values=c("black","grey80"), name="") +
  facet_wrap(~sp, scale="free", nrow=1) +
  theme_bw() + theme(legend.position="top")
r3
ggsave(r3, filename="Figures_Nmix/nmix_comparison_pdecay_all.tiff", height=6, width=8, units="in", dpi=300, compression="lzw")

r4 <- ggplot(data=ratio3) +
  geom_point(aes(x=N, y=scaledN, group=sp), pch=1, col="grey80") +
  geom_point(aes(x=N, y=scaled, group=sp, col=converged), size=3) +
  stat_smooth(aes(x=N, y=scaledN), method="lm", formula= y~x, col="grey80") + 
  stat_smooth(aes(x=N, y=scaled), method="lm", formula= y~x) + 
  scale_color_manual(values=c("black","grey80"), name="") +
  facet_wrap(~sp, scale="free", nrow=1) +
  theme_bw() + theme(legend.position="top")
r4
ggsave(r4, filename="Figures_Nmix/nmix_comparisonscaled_pdecay_all.tiff", height=6, width=8, units="in", dpi=300, compression="lzw")

#############
#4. look at output of all models
#############
require(tidyr)

# o.nmix <- data.frame(fullfilename=dir("Data_Model_Output/Output_Nmix/UI_base_t_2020-07-07/", pattern=".txt", full.names = T))
# o.nmix <- rbind(o.nmix, data.frame(fullfilename=dir("Data_Model_Output/Output_Nmix/UI_psite_t_2020-07-07/", pattern=".txt", full.names = T)))
# o.nmix <- rbind(o.nmix, data.frame(fullfilename=dir("Data_Model_Output/Output_Nmix/UI_pdecay_t_2020-07-07/", pattern=".txt", full.names = T)))
o.nmix <- data.frame(fullfilename=dir("Data_Model_Output/Output_Nmix/alldays_base/", pattern=".txt", full.names = T))
o.nmix <- rbind(o.nmix, data.frame(fullfilename=dir("Data_Model_Output/Output_Nmix/alldays_psite/", pattern=".txt", full.names = T)))
o.nmix <- rbind(o.nmix, data.frame(fullfilename=dir("Data_Model_Output/Output_Nmix/alldays_pdecay/", pattern=".txt", full.names = T)))

o.nmix <- separate(o.nmix, col="fullfilename", into=c("Data_Model_Output","outfolder","model","filename"), sep="/", remove=F)
o.nmix <- separate(o.nmix, col="filename", into=c("outjagsUI","mod","params","sp","grid","t","date"), sep="_", remove=F)
o.nmix$date <- gsub(o.nmix$date, pattern=".txt", replacement="")
o.nmix$ctime <- gsub(o.nmix$t, pattern="t", replacement="")

out.all <- data.frame()
for(o in o.nmix$fullfilename)
{
  print(o)
  
  out <- read.table(as.character(o), sep=",")
  
  out$mod <- o.nmix[o.nmix$fullfilename==o,]$mod
  out$date <- o.nmix[o.nmix$fullfilename==o,]$date
  out$sp <- o.nmix[o.nmix$fullfilename==o,]$sp
  out$grid <- o.nmix[o.nmix$fullfilename==o,]$grid
  out$ctime <- o.nmix[o.nmix$fullfilename==o,]$ctime
  out$param <- row.names(out)
  
  out.all <- rbind(out.all, out)
}

out.all$converged <- "not converged"
out.all[!is.na(out.all$Rhat) & out.all$Rhat < 1.1,]$converged <- "converged"

# write.table(out.all, file="Data_Final/output_Nmix_2021-06-28.txt", sep=",")
write.table(out.all, file="Output_cam/output_Nmix_2023-09-26.txt", sep=",")
