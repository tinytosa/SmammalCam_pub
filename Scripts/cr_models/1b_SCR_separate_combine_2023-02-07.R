#combine information for all SCR-B separate models for all grids and all species


require(tidyr)

SCRB.files <- data.frame(f = dir("C:/Users/tosam/Documents/0_OSU/Dissertation/Manuscripts/SmallMammalCameras/SmammalCam_pub/Output_SCR/", pattern="SCR_B", full.names=T))

SCRB.files$files <- gsub(SCRB.files$f, pattern="C:/Users/tosam/Documents/0_OSU/Dissertation/Manuscripts/SmallMammalCameras/SmammalCam_pub/", replacement="")

SCRB.files <- separate(SCRB.files, col="files", into=c("Output","filename"), sep="/", remove=F)
SCRB.files <- separate(SCRB.files, col="filename", into=c("sp","SMgrid","model","B","date"), sep="_", remove=F)
SCRB.files$grid <- gsub(SCRB.files$SMgrid, pattern="SM-0", replacement="")


scr.ests <- data.frame()
for(r in 1:nrow(SCRB.files))
{
  f <- as.character(SCRB.files$f[r])
  print(f)
  #read in entire file
  info <- read.delim(f, skipNul=T) #JAGS output? #S.Area <- nrow(data$ssDF[[1]])*data$res should be S.Area <- nrow(data$ssDF[[1]])*data$res^2
  
  first <- which(grepl(info$X.mod, pattern="n.eff$") == T) #gives row with header
  last <- which(grepl(info$X.mod, pattern="^deviance") == T) #gives last row of table
  
  info2 <- read.table(f, skipNul=T, skip=first, nrows=last-first) #read in just the parameter estimates
  info2$param <- row.names(info2)

  n.num <- which(info$X.mod=="$Nhat") #get number row for N
  n <- as.character(info[n.num+1,])
  if(grepl(n, pattern="var1"))
  {
    n <- as.character(info[n.num+2,])
  }
  n <- gsub(n, pattern="\\[1\\] ", replacement="")
  n <- trimws(n)
  
  d.num <- which(info$X.mod=="$Dhat") #get number row for D
  d <- as.character(info[d.num+1,])
  if(grepl(d, pattern="var1"))
  {
    d <- as.character(info[d.num+2,])
  }
  d <- gsub(d, pattern="\\[1\\] ", replacement="")
  d <- trimws(d)
  
  info2$sp <- SCRB.files$sp[r]
  info2$grid <- SCRB.files$grid[r]
  info2$mode <- NA
  info2[info2$param == "N",]$mode <- n
  info2[info2$param == "D",]$mode <- d
  
  scr.ests <- rbind(scr.ests, info2)
}

write.table(scr.ests, file="Output/SCR-B_separate.txt", sep=",")
