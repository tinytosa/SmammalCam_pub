#small mammal camera data
#cleaned up from SMcams_ProportionCamerasDetected_2019-10-28.R
#Charlotte Eriksson's method

#######################
#load packages
require(plyr) #for ddply
require(lubridate) #for dates
require(reshape2) #for melt

require(ggplot2)
require(ggpubr)

#
path.local <- getwd()

###############
#get sunrise and sunset times
#from https://www.esrl.noaa.gov/gmd/grad/solcalc/
#44.2129, -122.2552 (HJA headquarters)

sunrise <- read.csv("Data_Raw/2017_sunrise_noaa.csv")
sunrise <- melt(sunrise, id.vars=c("Day"))
sunrise <- sunrise[!sunrise$value=="",] #remove days that don't exist
names(sunrise) <- c("date","month","sunrise")
sunrise$sunrisedt <- as.POSIXct(paste(sunrise$month, sunrise$date, "2017", sunrise$sunrise , sep=" "), format="%b %d %Y %H:%M", tz="GMT")

sunset <- read.csv("Data_Raw/2017_sunset_noaa.csv")
sunset <- melt(sunset, id.vars=c("Day"))
sunset <- sunset[!sunset$value == "",] #remove days that don't exist
names(sunset) <- c("date","month","sunset")
sunset$sunsetd <- as.POSIXct(paste(sunset$month, sunset$date, "2017", sep=" "), format="%b %d %Y", tz="GMT")
sunset$sunsetdt <- as.POSIXct(paste(sunset$month, sunset$date, "2017", sunset$sunset , sep=" "), format="%b %d %Y %H:%M", tz="GMT")

#########
#trap information
traps <- read.table("Data_Raw/traplocs.txt", sep=",", header=T) #164 locs (100 sherman + 64 tomahawk)
pema.grid <- traps[traps$PEMAgrid == 1,] #25 locs
tato.grid <- traps[traps$ttype == "tomahawk",] #64 locs

#keywords to species and trap type
ks <- data.frame(Keywords=c("FLYING SQUIRREL","DEER MOUSE","TOWNSENDS CHIPMUNK","VOLE"), sp=c("GLOR","PEMA","TATO","MYO"), ttype=c("tomahawk","sherman","tomahawk","sherman"))

#########
#load data in from here
#########
cam.data <- read.table("Data_Raw/SMcamera_tags_daynum_2020-05-04.txt", sep=",")
cam.data$DT <- as.POSIXct(cam.data$DateTimeOriginal, format="%Y-%m-%d %H:%M:%S", tz="GMT")
cam.data <- cam.data[cam.data$grid != 0,] #remove practice grid

#determine which cameras were active before removing detections
active.cams <- unique(cam.data[,c("Directory","grid","Loc")]) #number of active cameras
active.cams <- merge(active.cams, traps[,c("RC","ttype","PEMAgrid")], all.x=T, by.x="Loc", by.y="RC")
active.sherman <- active.cams[active.cams$PEMAgrid == 1,]
# active.sherman <- data.frame(table(active.cams[active.cams$PEMAgrid == 1,"grid"])) #can't summarize by grid because need to calculate standard deviation/error
active.tomahawk <- active.cams[active.cams$ttype == "tomahawk",]
# active.tomahawk <- data.frame(table(active.cams[active.cams$ttype == "tomahawk","grid"]))

#calculate number of days each camera was active by taking max (day) num for each camera image
max.active.cams <- ddply(cam.data, .(Directory, grid, Loc, Rlet, Cnum), summarize, maxday=max(daynum))
max.active.cams[max.active.cams$maxday < 5,]

#remove detections from sp not in analysis
cam.data <- merge(cam.data, ks, by="Keywords", all.x=T) #attach sp info to each line
cam.data <- cam.data[!cam.data$Keywords %in% "VOLE",] #remove voles
cam.data <- cam.data[(cam.data$sp %in% c("GLOR","TATO") & cam.data$Loc %in% tato.grid$RC) | (cam.data$sp == "PEMA" & cam.data$Loc %in% pema.grid$RC),] #remove sherman from GLSA and TATO, remove tomahawk from PEMA

#######
#summarize cam.data by first day of detection per species
day.data <- ddply(cam.data, .(grid, Directory, Loc, Rlet, Cnum, Keywords, sp, num), summarize, num.photos=length(Keywords)) #summarize number of detections per (day) num
day.data <- ddply(day.data, .(grid, Loc, Keywords, sp), summarize, first=min(num), totaldets=sum(num.photos)) #summarize by first day of detection
day.data$day <- as.integer(day.data$first)
day.data$grid <- factor(day.data$grid)
day.data$Directory <- paste(day.data$grid, day.data$Loc, sep="")

day.data <- merge(day.data, traps[,c("RC","posx","posy")], by.x="Loc", by.y="RC", all.x=T)

p.pema <- ggplot(data=day.data[day.data$sp == "PEMA",]) + geom_point(aes(x=posx, y=posy, size=totaldets)) + theme_bw(base_size=20) + facet_wrap(~grid, nrow=3) + theme(panel.grid = element_blank()) + coord_equal()
p.tato <- ggplot(data=day.data[day.data$sp == "TATO",]) + geom_point(aes(x=posx, y=posy, size=totaldets)) + theme_bw(base_size=20) + facet_wrap(~grid, nrow=3) + theme(panel.grid = element_blank()) + coord_equal()
p.glsa <- ggplot(data=day.data[day.data$sp == "GLOR",]) + geom_point(aes(x=posx, y=posy, size=totaldets)) + theme_bw(base_size=20) + facet_wrap(~grid, nrow=3) + theme(panel.grid = element_blank()) + coord_equal()

# ggsave(ggarrange(p.pema, p.tato, p.glsa, nrow=1, align="h"), file="Figures/cam_totaldetectons.tiff", height=8, width=24, units="in", dpi=400, compression="lzw")

#####
#no subsetting
#number of cameras species was detected on if camera number is not important
#don't do this
#####
# table(day.data[,c("grid","Keywords")])
#grid DEER MOUSE FLYING SQUIRREL TOWNSENDS CHIPMUNK
#   0         25               8                 25 #practice grid
#   1         71              66                 76
#   2         69              46                 78
#   3         75              73                 77
#   4         75              61                 68
#   5         76              59                 77
#   6         76              49                 67
#   7         77              53                 80
#   8         75              64                 78

# after removing sherman grid from TATO and GLSA, remove tomahawk grid from PEMA
# grid DEER MOUSE FLYING SQUIRREL TOWNSENDS CHIPMUNK
#    1         18              56                 62
#    2         20              39                 62
#    3         22              58                 62
#    4         24              48                 52
#    5         23              45                 62
#    6         25              38                 54
#    7         24              37                 64
#    8         24              50                 62

# detections <- data.frame(table(day.data[,c("grid","sp")])) #change this endpoint summary data

#what if no detection of species at a camera?
#need to check if camera had photos first
#then fill in cameras with no detection of species

# detections$operational <- c(active.tomahawk$Freq, active.sherman$Freq, active.tomahawk$Freq)
# detections$prop.cam.detected <- detections$Freq/detections$operational
# 
# detections %>% group_by(sp) %>% dplyr::summarize(mean = mean(prop.cam.detected))
#   sp     mean
# 1 GLOR  0.733
# 2 PEMA  0.913
# 3 TATO  0.949

#############################################
#consolidate detections based on consolidation window

avg.det.vals <- data.frame()
for(t in c(0, 15, 60, 1440)) #t in minutes
{
  print(t)
  for(d in 1:9) 
  {
    print(d)
    data <- cam.data[cam.data$num <= d,] #up to day number d
    
    for(sp in c("PEMA","GLOR","TATO"))
    {
      ttype <- ks[ks$sp == sp,]$ttype
      if(ttype == "sherman")
      {
        active <- active.sherman
      }
      
      if(ttype == "tomahawk")
      {
        active <- active.tomahawk
      }
      
      #order the photos by date and time
      data.sp <- data[data$sp == sp,]
      data.sp <- data.sp[order(data.sp$DT),]
      data.sp <- data.sp[order(data.sp$Directory),]
      #consolidate photos within t minutes
      data.sp$timediff <- c(0,round((data.sp$DT[-1] - data.sp$DT[-nrow(data.sp)])/60, digits=0)) #need to do this for each species
      n <- 1
      data.sp$consolidate <- sapply(data.sp$timediff, function(x) {if(x > t | x < 0) #edited this 7/6/2020, make sure to run again
      {n <<- n+1} 
        n})
      
      #consolidate if keywords and consolidate numbers are the same
      data.sp <- ddply(data.sp, .(Directory, grid, Loc, Rlet, Cnum, sp, num, consolidate, nocturnal), summarise, numphoto=length(DT), startDT=min(DT))
      
      ######
      avg.det <- ddply(data.sp, .(Directory, grid, sp), nrow) #number of detections per grid
      # avg.det <- merge(avg.det, active, by.x="grid", by.y="Var1") #merge with number of number of active cameras per grid, can't do this if trying to calculate error bars
      # avg.det$prop <- avg.det$V1/avg.det$Freq #calculate average detections per camera per grid, can't do this anymore since each camera is represented in the data
      # avg.det$t <- t
      # avg.det$d <- d
      
      avg.det <- merge(avg.det, active, by=c("Directory","grid"), all.y=T) #can check against Figures/cam-totaldetecitions.tiff. missing points also reflect non-operational cameras
      avg.det[is.na(avg.det$V1),]$sp <- sp #fix values for cameras with no detections of species
      avg.det[is.na(avg.det$V1),]$V1 <- 0
      avg.det.sp <- ddply(avg.det, .(grid, sp), summarize, mean=mean(V1), sd=sd(V1), ncams=length(V1))
      avg.det.sp$t <- t
      avg.det.sp$d <- d
      
      avg.det.vals <- rbind(avg.det.vals, avg.det.sp)
    }
  }
}

# write.table(avg.det.vals, file="Output/output_avgdets_2023-03-10.txt", sep=",")
