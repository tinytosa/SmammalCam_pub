#time to event models
#stations should not be baited#

#install.packages("remotes")
#remotes::install_github("annam21/spaceNtime")
#Sys.unsetenv("GITHUB_PAT")
#Sys.getenv("GITHUB_PAT")
#devtools::install_github("annam21/spaceNtime", force = T, build_vignettes = T)
#browseVignettes("spaceNtime")
require(spaceNtime)

require(plyr) #for ddply
require(tidyr) #for separate
require(lubridate) #for dates
require(reshape2) #for melt
require(ggplot2)
require(ggpubr)
require(xlsx)

path.local <- "C:/Users/tosam/Documents/0_OSU/Dissertation/Data/SmallMammalCams/"
#path.local <- "/raid1/home/fw/tosam/SmallMammalCams/"
setwd(path.local)

########
#import data
########
smdataall <- read.table("SMcamera_tags_2020-03-25.txt", sep=",") #all small mammal data
table(smdata[,c("Keywords","grid")])

smdata <- read.table("SMcamera_tags_daynum_2020-05-04.txt", sep=",") #only inlcudes data for deer mice, flying squirrels, townsends chipmunk, and voles
smdata$DT <- as.POSIXct(smdata$DT, format="%Y-%m-%d %H:%M:%S", tz="GMT")
smdata$D <- as.POSIXct(smdata$D, format="%Y-%m-%d", tz="GMT")

smdata <- smdata[smdata$num > 3 & smdata$num < 9,] #limit data to days when no bait present

traplocs <- read.table("traplocs.txt", sep=",", as.is = T) #as.is imports all columns with characters as character class

traps <- unique(smdata[,c("Directory","grid")])
traps$Rlet <- substr(traps$Directory, 2, 2)
traps.t <- traps[!traps$Rlet %in% c("P","R","T","V","X"),] #tomahawk traps

########
#cam info for deploy date times
########
#don't use these for TATO since start from "unbaited" days
# caminfo <- read.xlsx("checkphotos_2019-12-17.xlsx", sheetName="SmammalCamInfo_2019-12-21")
# caminfo <- caminfo[!caminfo$TrapName %in% c("4P3","5R3","6R1"),] #mismatches in date set and photo date times
# caminfo <- caminfo[order(caminfo$TrapName),]
# caminfo <- separate(caminfo, col=DateTimeBaited, into=c("DateBaited"), remove=F, sep=" ")
# caminfo$DateBaited <- as.POSIXct(caminfo$DateBaited, format="%Y-%m-%d", tz="GMT")

#use these for TATO, PEMA, GLSA
sunset <- read.csv("2017_sunset_noaa.csv")
sunset <- melt(sunset, id.vars=c("Day"))
sunset <- sunset[!sunset$value == "",] #remove days that don't exist
names(sunset) <- c("date","month","sunset")
sunset$sunsetd <- as.POSIXct(paste(sunset$month, sunset$date, "2017", sep=" "), format="%b %d %Y", tz="GMT")
sunset$sunsetdt <- as.POSIXct(paste(sunset$month, sunset$date, "2017", sunset$sunset , sep=" "), format="%b %d %Y %H:%M", tz="GMT")

sunrise <- read.csv("2017_sunrise_noaa.csv")
sunrise <- melt(sunrise, id.vars=c("Day"))
sunrise <- sunrise[!sunrise$value=="",] #remove days that don't exist
names(sunrise) <- c("date","month","sunrise")
sunrise$sunrised <- as.POSIXct(paste(sunrise$month, sunrise$date, "2017", sep=" "), format="%b %d %Y", tz="GMT")
sunrise$sunrisedt <- as.POSIXct(paste(sunrise$month, sunrise$date, "2017", sunrise$sunrise , sep=" "), format="%b %d %Y %H:%M", tz="GMT")

#set sunset times to DateBaited for nocturnal species
#caminfo <- merge(caminfo, sunset[,c("sunsetd","sunsetdt")], all.x=T, by.x="DateBaited", by.y="sunsetd")

########
#all camera active days
smdataall <- separate(smdataall, col="DateTimeOriginal", into=c("Date"), remove=F, sep=" ")
caminfo.days <- ddply(smdataall, .(Directory, grid, Loc, Rlet, Cnum, Date), nrow)
caminfo.days$Date <- as.POSIXct(caminfo.days$Date, format="%Y-%m-%d", tz="GMT")
caminfo.days <- merge(caminfo.days, sunset[,c("sunsetd","sunsetdt")], by.x="Date", by.y="sunsetd", all.x=T)

caminfo.days <- merge(caminfo.days, sunrise[,c("sunrised","sunrisedt")], by.x="Date", by.y="sunrised", all.x=T)
caminfo.days$NextDate <- caminfo.days$Date + days(1)
caminfo.days <- merge(caminfo.days, sunrise[,c("sunrised","sunrisedt")], by.x="NextDate", by.y="sunrised", all.x=T, suffixes=c(".diu",".noct"))

#area
r <- 5 # camera maximum distance in meters
theta <- 44 # camera viewshed in degrees
a <- pi * r^2 * theta/360 # Area of a single camera in square meters

#parameters
params <- data.frame(sp=c("DEER MOUSE","TOWNSENDS CHIPMUNK","FLYING SQUIRREL"),
                     t=c(15,15,15),
                     m=c(24.09, 86.58, 87.22), #derived from mmdm from allgrids SCR models
                     sa=rep(280*280, 3),
                     active=c("night","day","night"))

results <- data.frame()
for(sp in c("DEER MOUSE","TOWNSENDS CHIPMUNK","FLYING SQUIRREL"))
{
  print(sp)
  
  #t <- params[params$sp == sp,]$t # duration of each sampling occasion for STE
  m <- params[params$sp == sp,]$m # animal speed in m per hour for TTE
  sa <- params[params$sp == sp,]$sa # study area size
  print(paste(m, sa))
  
  for(g in c(1:8))
    {
    print(g)
    sp.g.df <- smdata[smdata$grid == g & smdata$Keywords == sp,]
    sp.g.df <- sp.g.df[order(sp.g.df$DT),]
    sp.g.df <- sp.g.df[order(sp.g.df$Directory),]
    
    deploydates <- ddply(sp.g.df, .(grid), summarise, start=min(D), end=max(D))
    
    ####### wrong, don't use?
    #cut into occasions, don't consolidate, photos within t minutes
    # sp.g.df$timediff <- c(0,round((sp.g.df$DT[-1] - sp.g.df$DT[-nrow(sp.g.df)])/60, digits=0)) #need to do this for each species
    # 
    # n <- 1
    # sp.g.df$consolidate <- sapply(sp.g.df$timediff, function(x) {if(x > t | x < 0)
    # {n <<- n+1} 
    #   n})
    # 
    # #consolidate if keywords and consolidate numbers are the same
    # con.df <- ddply(sp.g.df, .(Directory, grid, Loc, Rlet, Cnum, Keywords, num, consolidate, nocturnal), summarise, numphoto=length(DT), startDT=min(DT))
    #######
    
    #time to first detection (occasion is whole sampling period)
    #con.df <- ddply(sp.g.df, .(Directory, grid, Loc, Rlet, Cnum, Keywords), summarise, startDT=min(DT))
    
    #or first detection per occasion by num (day or night num)
    #con.df <- ddply(sp.g.df, .(Directory, grid, Loc, Rlet, Cnum, Keywords, num), summarise, startDT=min(DT))
    
    #or use all detections
    con.df <- sp.g.df
    con.df$startDT <- con.df$DT
    
    con.df <- con.df[con.df$Directory %in% unique(caminfo.days$Directory),] #make sure Directory is in caminfo
    
    df <- data.frame(cam = con.df$Directory,
                     datetime = con.df$startDT,
                     count = rep(1, nrow(con.df))) #for just present/absent
    #count = con.df$numphoto)
    
    if(params[params$sp == sp,]$active == "night")
    {
      #one date time for each camera for nocturnal species
      deploy <- data.frame(cam = unique(caminfo.days[caminfo.days$grid == g,]$Directory),
                           start = rep(sunset[sunset$sunsetd == deploydates$start,]$sunsetdt, 
                                       length(unique(caminfo.days[caminfo.days$grid == g,]$Directory))), #for nocturnal
                           end = rep(sunrise[sunrise$sunrised == (deploydates$end + days(1)),]$sunrisedt, 
                                     length(unique(caminfo.days[caminfo.days$grid == g,]$Directory))), #for nocturnal
                           area = rep(a, length(unique(caminfo.days[caminfo.days$grid == g,]$Directory))))
    }
    else
    {
      #one date time for each camera for diurnal species
      deploy <- data.frame(cam = unique(caminfo.days[caminfo.days$grid == g,]$Directory),
                           start = rep(sunrise[sunrise$sunrised == deploydates$start,]$sunrisedt, 
                                       length(unique(caminfo.days[caminfo.days$grid == g,]$Directory))), #for diurnal
                           end = rep(sunset[sunset$sunsetd == deploydates$end,]$sunsetdt, 
                                     length(unique(caminfo.days[caminfo.days$grid == g,]$Directory))), #for diurnal
                           area = rep(a, length(unique(caminfo.days[caminfo.days$grid == g,]$Directory))))
    }

    #######
    #space-to-event for motion detection photos formatted like time lapse photos
    #######
    for(t in c(1, 5, 15, 30, 60))
      {
      study_dates <- c(min(deploy$start), max(deploy$end))
      occ <- build_occ(samp_freq = 3600, # seconds between the start of each sampling occasion = 1 hr
                       samp_length = 60*t, # 60*number of minutes: duration of each sampling occasion (seconds)
                       study_start = study_dates[1],
                       study_end = study_dates[2])
      if(params[params$sp == sp,]$active == "night")
        {
        occ <-  occ[hour(occ$start) > 17 | hour(occ$start) < 8,] #remove day time hours for nocturnal, hours 8-17
        }
      else
        {
        occ <-  occ[hour(occ$start) > 6 & hour(occ$start) < 19,] #remove night time hours for diurnal, hours before 7 and after 18
        }
      
      ste_eh <- ste_build_eh(df, deploy, occ)
      if(all(is.na(ste_eh$STE)))
        {
        ste <- data.frame(N=NA, SE=NA, LCI=NA, UCI=NA, model="ste")
        }
      else
        {
        ste <- ste_estN_fn(ste_eh, study_area = sa) # 280m*280m grid
        #names(ste) <- paste(names(ste), "ste", sep="_")
        ste$model <- "ste"
        }
      
      results <- rbind(results, data.frame(sp, g, t, m, a, sa, ste))
      }
      
      #######
      #time-to-event
      #######
      occ.tte <- build_occ(samp_freq = 3600 * 24, # 3600*number of hours
                           samp_length = 3600 * 16, # 3600*number of hours
                           study_start = study_dates[1],
                           study_end = study_dates[2])
      
      #mean amount of time (in seconds) that it takes for an animal to cross the average viewshed of a camera
      per <- tte_samp_per(deploy, lps = m/3600) # lps is animal speed in m/sec, could use average distance moved?
      
      tte_eh <- tte_build_eh(df, deploy, occ.tte,  per)
      tte <- tte_estN_fn(tte_eh, study_area = sa) # 280m*280m grid
      #names(tte) <- paste(names(tte), "tte", sep="_")
      tte$model <- "tte"
      
      results <- rbind(results, data.frame(sp, g, t=1440, m, a, sa, tte))
  }
}

write.table(results, file=paste("STE_TTE_model_results_", Sys.Date(), ".txt", sep=""), sep=",", row.names=F)

#count ~ pois(lambda)
#time ~ exp(lambda)

############
#post-processing ste and tte
############

require(ggplot2)
require(ggpubr)

sptospecies <- data.frame(sp=c("PEMA","TATO","GLSA"),
                          species=c("DEER MOUSE","TOWNSENDS CHIPMUNK", "FLYING SQUIRREL"))

###########
#import results from oSCR allgrids
###########
# oscr <- read.table("Output_oSCR_allgrids_2020-08-26.txt", sep=",")
# oscr$g <- substr(oscr$grid, 5, 5)
# oscr$sp.y <- oscr$sp
# oscr$scaledscrn <- oscr$scaledn
#oscr$model.x <- "oSCR all"
#fixed in SMcamera_postprocessing_SCRdata
# oscr$pixels <- oscr$n/oscr$estimate
# oscr$n_lwr <- oscr$lwr*oscr$pixels
# oscr$n_upr <- oscr$upr*oscr$pixels

###########
#import results from SCR random
#scaledscrn is incorrect
###########
oscr <- read.table("Output_SCR_randomeff_2020-08-12.txt", sep=",", header=T)
for(sp in unique(oscr$sp))
{
  print(sp)
  toscale <- oscr[oscr$sp == sp,]
  oscr[oscr$sp == sp,]$scaledscrn <- scale(toscale$nmode)
}
oscr$sp.y <- factor(oscr$sp, levels=c("PEMA","TATO","GLSA"), labels=c("PEMA","TATO","GLOR"))
oscr$g <- oscr$grid
oscr$mod <- oscr$model


eventmodel <- read.table("STE_TTE_model_results_2020-08-26.txt", sep=",", header=T)
eventmodel <- merge(eventmodel, sptospecies, all.x=T, by.x="sp", by.y="species")
eventmodel$model <- factor(toupper(eventmodel$model), levels=c("TTE","STE"))
eventmodel$sp.y <- factor(eventmodel$sp.y, levels=c("PEMA","TATO","GLSA"), labels=c("PEMA","TATO","GLOR"))
eventmodel$scaledn <- 0

for(sp in unique(eventmodel$sp.y))
{
  for(f in unique(eventmodel$t))
  {
    print(paste(sp, f))
    toscale <- eventmodel[eventmodel$sp.y == sp & eventmodel$t == f,]
    eventmodel[eventmodel$sp.y == sp & eventmodel$t == f,]$scaledn <- scale(toscale$N)
  }
}

plotdata <- merge(eventmodel, oscr[,!names(oscr) %in% c("sp","sp.1","model")], all.x=T, by.x=c("g", "sp.y"), by.y=c("g", "sp.y"))
#plotdata$sp.y <- factor(plotdata$sp.y, levels=c("PEMA","TATO","GLSA"), labels=c("PEMA","TATO","GLOR"))
plotdata$t <- factor(plotdata$t)
plotdata$grid <- paste("SM-0", plotdata$g, sep="")
plotdata$mod <- plotdata$model


#plot by scaled estimates
est_comp <- ggplot(data=plotdata, aes(x=nmode, y=N, col=t, shape=model)) + 
  geom_point(size=3) +
  #geom_errorbar(aes(ymin=LCI, ymax=UCI), width=0.1) +
  #geom_errorbar(aes(xmin=n_lwr, xmax=n_upr), width=0.1) +
  stat_smooth(aes(group=t), method="lm", formula=y~x, fill=NA) +
  #stat_smooth(aes(x=nmode, y=scaledscrn), method="lm", formula= y~x, col="grey80", fill=NA, lwd=1) +
  #xlab("oSCR_all estimate") +
  xlab("SCR random estimate (N)") +
  ylab("unmarked model estimate") +
  #geom_text(aes(label=grid)) +
  facet_wrap(~model + sp.y, scale="free", nrow=2) + theme_bw(base_size = 15)
est_comp
#ggsave(est_comp, filename=paste("Figures_TTE_STE/ModelEstimates_TTE_STE_vs_oSCR_", Sys.Date(), ".tiff", sep=""), 
ggsave(est_comp, filename=paste("Figures_TTE_STE/ModelEstimates_TTE_STE_vs_SCRrandom_", Sys.Date(), ".tiff", sep=""), 
       width=10, height=8, dpi=300, compression="lzw")

#plot by grid
est <- ggplot(data=plotdata, aes(x=g, col=t, shape=model)) +
  geom_point(data=oscr[!oscr$g ==9,], aes(y=n, x=grid), col="black") +
  #geom_errorbar(data=oscr[!oscr$g ==9,], aes(xmin=n_lwr, xmax=n_upr, y=grid), col="black", width=0.1) + #for oSCRall estimates
  geom_errorbar(data=oscr[!oscr$g ==9,], aes(ymin=n95low, ymax=n95high, y=grid), col="black", width=0.1) + #for SCRrandom estimates
  
  geom_point(aes(y=N), size=3, position=position_dodge(width=0.75)) +
  geom_errorbar(aes(ymin=LCI, ymax=UCI), width=0.1, position=position_dodge(width=0.75)) +
  #scale_x_log10() +
  coord_flip() +
  facet_wrap(~sp.y, scale="free") + theme_bw(base_size=15)
est
ggsave(est, filename=paste("Figures_TTE_STE/ModelEstimates_TTE_STE_", Sys.Date(), ".tiff", sep=""), 
       width=10, height=6, dpi=300, compression="lzw")

#######
#scale and run linear regression
#######

#scale values by species and model type and compare to SCR_random results
plotdata$scaledn <- 0
slopes <- data.frame()
for(sp in unique(plotdata$sp.y))
{
  for(m in c("STE"))
  {
    for(t in c(1,5,15,30,60))
    {
      print(sp); print(m); print(t)
      
      plotdata[plotdata$sp.y == sp & plotdata$mod == m & plotdata$t == t,]$scaledn <- scale(plotdata[plotdata$sp.y == sp  & plotdata$mod == m & plotdata$t == t,]$N)
      
      #build in linear regressions here
      lm.scaled <- lm(data=plotdata[plotdata$sp.y == sp & plotdata$mod == m & plotdata$t == t,], scaledn ~ scaledscrn)
      rmse <- sqrt(mean(lm.scaled$residuals^2))
      #sink(paste("Output_TTE_STE/scaled_lm_", sp, "_",m,"_t", t, "_oscrall.txt", sep=""))
      sink(paste("Output_TTE_STE/scaled_lm_", sp, "_",m,"_t", t, "_scrrandom.txt", sep=""))
      print(summary(lm.scaled))
      sink()
      slopes <- rbind(slopes, data.frame(sp=sp, mod=m, t=t, slope=lm.scaled$coefficients[2], sd=coef(summary(lm.scaled))[,"Std. Error"][2], rmse=rmse))
      
    }
  }
  for(m in c("TTE"))
  {
    print(sp); print(m)
    
    plotdata[plotdata$sp.y == sp & plotdata$mod == m,]$scaledn <- scale(plotdata[plotdata$sp.y == sp  & plotdata$mod == m,]$N)
    
    #build in linear regressions here
    lm.scaled <- lm(data=plotdata[plotdata$sp.y == sp & plotdata$mod == m,], scaledn ~ scaledscrn)
    rmse <- sqrt(mean(lm.scaled$residuals^2))
    #sink(paste("Output_TTE_STE/scaled_lm_", sp, "_",m,"_oscrall.txt", sep=""))
    sink(paste("Output_TTE_STE/scaled_lm_", sp, "_",m,"_scrrandom.txt", sep=""))
    print(summary(lm.scaled))
    sink()
    slopes <- rbind(slopes, data.frame(sp=sp, mod=m, t=1440, slope=lm.scaled$coefficients[2], sd=coef(summary(lm.scaled))[,"Std. Error"][2], rmse=rmse))
    
    
  }
}

#plot TTE and STE scaled estimates vs. SCR_random/oSCR_all estimates
#no oSCR data for deer mouse grid 6 and STE data for GLSA t=1 for grid 2 and 5
p1 <- ggplot(data=plotdata, aes(x=nmode, y=scaledn, col=t, shape=mod)) +
  geom_point(size=3) +
  geom_smooth(aes(x=nmode, y=scaledscrn), method="lm", col="grey70", lwd=0.25) +
  geom_smooth(aes(col=t, lty=mod), method="lm", lwd=1, fill=NA) +
  #geom_text(aes(label=grid), nudge_y = 0.1) +
  # scale_color_manual(breaks=c("minimum known alive","Huggins Null", "oSCR0","oSCR all",
  #                             "SCR-B trapslarge","SCR-B allgrids","SCR-B randomeff", "SC0"),
  #                    labels=c("MNKA","Huggins","oSCR0","oSCR all",
  #                             "SCR0","SCR all","SCR random","SC0"),
  #                    values=c("black","grey70","#83B6C3","steelblue",
  #                             "#8C2155","#B9314F","#FCBFB7","#AAC0AA")) +
  # scale_shape_manual(breaks=c("minimum known alive","Huggins Null", "oSCR0","oSCR all",
  #                             "SCR-B trapslarge","SCR-B allgrids","SCR-B randomeff", "SC0"),
  #                    labels=c("MNKA","Huggins","oSCR0","oSCR all",
  #                             "SCR0","SCR all","SCR random","SC0"),
  #                    values=c(19,15,17,17,18,18,18,1)) +
  # scale_linetype(breaks=c("minimum known alive","Huggins Null", "oSCR0","oSCR all",
  #                         "SCR-B trapslarge","SCR-B allgrids","SCR-B randomeff", "SC0"),
  #                labels=c("MNKA","Huggins","oSCR0","oSCR all",
  #                         "SCR0","SCR all","SCR random","SC0")) +
  ylab("Scaled model estimate") + #xlab("oSCR all estimate (N)") + 
  xlab("SCR random all estimate (N)") +
  facet_wrap(~ mod + sp.y, scales="free") + theme_bw(base_size=15)
p1
ggsave(p1, filename=paste("Figures_TTE_STE/ModelEstimates_vs_SCRrandom_", Sys.Date(), ".tiff", sep=""), width=12, height=8, units="in", dpi=300, compression="lzw")

#plot scaled estimate vs. SCR_B estimate slopes
#slopes$mod <- factor(slopes$mod, levels=c("minimum known alive","Huggins Null", "oSCR0","oSCR all",
#                                          "SCR-B trapslarge","SCR-B allgrids","SCR-B randomeff", "SC0"))
slopes$sp <- factor(slopes$sp, levels=c("PEMA","TATO","GLOR"))
slopes$t <- factor(slopes$t)
slopes$mod <- factor(slopes$mod, levels=c("TTE","STE"))
pos <- position_dodge(width=0.5)

p2 <- ggplot(data=slopes, aes(x=sp, y=slope, col=t, shape=mod)) +
  geom_hline(aes(yintercept=1), lty="dashed", col="grey30", lwd=1) +
  geom_hline(aes(yintercept=0), lty="dashed", col="grey30", lwd=1) +
  
  geom_point(size=3, position=pos) +
  geom_errorbar(aes(ymin=slope-sd, ymax=slope+sd), position=pos, width=0.1) +
  # scale_color_manual(breaks=c("minimum known alive","Huggins Null", "oSCR0","oSCR all",
  #                             "SCR-B trapslarge","SCR-B allgrids","SCR-B randomeff", "SC0"),
  #                    labels=c("MNKA","Huggins","oSCR0","oSCR all",
  #                             "SCR0","SCR all","SCR random","SC0"),
  #                    values=c("black","grey70","#83B6C3","steelblue",
  #                             "#8C2155","#B9314F","#FCBFB7","#AAC0AA")) +
  ylab("slope") +
  #scale_y_continuous(breaks=seq(-0.2, 1.2, by=0.2)) +
  #ylim(c(-0.2,1.2)) +
  theme_bw(base_size=15)
p2
#ggsave(p2, filename="Figures_TTE_STE/ModelEstimates_slopes_vs_oSCRall_2020-08-27.tiff", height=6, width=8, dpi=300, units="in", compression="lzw")
ggsave(p2, filename="Figures_TTE_STE/ModelEstimates_slopes_vs_SCRrandom_2021-03-30.tiff", height=6, width=8, dpi=300, units="in", compression="lzw")

#plot rmse
p3 <- ggplot(data=slopes, aes(x=sp, y=rmse, col=t, shape=mod)) +
  geom_point(size=3, position=pos) +
  # scale_color_manual(breaks=c("minimum known alive","Huggins Null", "oSCR0","oSCR all",
  #                             "SCR-B trapslarge","SCR-B allgrids","SCR-B randomeff", "SC0"),
  #                    labels=c("MNKA","Huggins","oSCR0","oSCR all",
  #                             "SCR0","SCR all","SCR random","SC0"),
  #                    values=c("black","grey70","#83B6C3","steelblue",
  #                             "#8C2155","#B9314F","#FCBFB7","#AAC0AA")) +
  ylim(c(0,1.0)) +
  ylab("RMSE") +
  theme_bw(base_size=15)
p3

#ggsave(ggarrange(p2, p3, ncol=1, common.legend=T, align="v", heights=c(2,1)), filename="Figures_TTE_STE/ModelEstimates_slopes_RMSE_vs_oSCR_2020-08-27.tiff", height=6, width=8, dpi=300, units="in", compression="lzw")
ggsave(ggarrange(p2, p3, ncol=1, common.legend=T, align="v", heights=c(2,1)), filename="Figures_TTE_STE/ModelEstimates_slopes_RMSE_vs_SCRrandom_2021-03-30.tiff", height=6, width=8, dpi=300, units="in", compression="lzw")
