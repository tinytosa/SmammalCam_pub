#combine estimates from all unmarked models

######
#load packages
#####
require(tidyr)
require(ggplot2)

#conversion table
sp.convert <- data.frame(sp=c("chipmunk","flying\nsquirrel","flying\nsquirrel","mouse"), acro=c("TATO","GLSA","GLOR","PEMA"))
# sp.convert$sp <- factor(sp.convert$sp, levels=c("mouse","chipmunk","flying\nsquirrel"))

###################
avgdets <- read.table("Output_cam/output_avgdets_2023-03-10.txt", sep=",")
avgdets <- avgdets[avgdets$grid != 0,] #remove practice grid if in data
table(avgdets$sp)

avgdets$t_d <- paste("avgdets_", "t", avgdets$t, "_d", avgdets$d, sep="")

avgdets$est <- avgdets$mean
avgdets$se <- avgdets$sd/sqrt(avgdets$ncams)
avgdets$D2.5 <- avgdets$mean - 1.96*avgdets$sd/sqrt(avgdets$ncams)
avgdets$D97.5 <- avgdets$mean + 1.96*avgdets$sd/sqrt(avgdets$ncams)
avgdets$model <- avgdets$t_d

avgdets$meansd <- paste(sprintf('%.2f',avgdets$mean), "\u00B1", sprintf('%.2f', avgdets$sd))
# avgdets$meansd <- paste(sprintf('%.2f', avgdets$est), " (", sprintf('%.2f', avgdets$D2.5), " - ", sprintf('%.2f', avgdets$D97.5), ")", sep="")

avgdets.wide <- pivot_wider(avgdets, id_cols=c('grid','sp'), names_from='t_d', values_from = 'meansd')
avgdets.wide[avgdets.wide$sp == "GLOR",]$sp <- "GLSA"

###################
#mmdm values are from TATO_cr_date.txt files in Output_SCR folder, SCR random data formatting
info <- data.frame(sp=c("GLSA","TATO","PEMA"), mmdm=c(84.11, 79.48, 24.09), length=c(280,280,90)) #length = length of each side of grid, 7*40, 9*10
a <- function(buffer, length) {pi*buffer^2 + buffer*length*4 + length*length} # area circle with r=3*mmdm + buffer*length of grid*4 + lengthofgrid*lengthofgrid

info$mmdm.a <- a(info$mmdm, info$length)/10000 #area in ha for 1 mmdm buffer around grid
info$halfmmdm.a <- a(info$mmdm/2, info$length)/10000 #area in ha for 0.5 mmdm buffer around grid
info$n.occ.ideal <- c(8,8,4)

info
#     sp  mmdm length    mmdm.a halfmmdm.a
# 1 GLSA 84.11    280 19.482837  13.105789
# 2 TATO 79.48    280 18.726326  12.787022
# 3 PEMA 24.09     90  1.859555   1.289199

#Nmixture model estimates
nmix <- read.table("Output_cam//output_Nmix_2023-09-26.txt", sep=",")
nmix <- nmix[nmix$param == "lambda" & nmix$grid > 0 & nmix$sp != "myo",]

nmix$sp <- toupper(nmix$sp)

#turn values into factors
# nmix <- merge(sp.convert, nmix, by.y="sp", by.x="acro", all.y=T) #factor species names
nmix$t <- factor(nmix$ctime, levels=c(0,15,60,1440)) #factor delta 
nmix$mod <- factor(nmix$mod, levels=c("base","psite","pdecay"), labels=c("base","pstation","pdecay")) #factor model

nmix <- merge(nmix, info[,c("sp","mmdm.a","halfmmdm.a")], by.x="sp", by.y="sp", all.x=T)
nmix$mean.density <- nmix$mean/nmix$halfmmdm.a
nmix$mean.densitysd <- paste(sprintf('%.2f', nmix$mean.density), " (", sprintf('%.2f', nmix$X2.5./nmix$halfmmdm.a), " - ", sprintf('%.2f', nmix$X97.5./nmix$halfmmdm.a), ")", sep="")
nmix$mod_t <- paste("nmix_", nmix$mod, "_t", nmix$t, sep="")

nmix$est <- nmix$mean.density
nmix$D2.5 <- nmix$X2.5./nmix$halfmmdm.a
nmix$D97.5 <- nmix$X97.5./nmix$halfmmdm.a
nmix$model <- nmix$mod_t

nmix.wide <- pivot_wider(nmix, id_cols=c('sp','grid'), names_from='mod_t', values_from='mean.densitysd')

###########################
tte <- read.table("Output_cam/output_STE_TTE_model_results_2020-08-26.txt", sep=",", header=T)
#convert to density
tte$D <- tte$N/(tte$sa/10000)
tte$X2.5. <- tte$LCI/(tte$sa/10000)
tte$X97.5. <- tte$UCI/(tte$sa/10000)

tte$d.ci <- paste(round(tte$D, digits=2), " (", sprintf('%.2f',tte$X2.5.), " - ", sprintf('%.2f', tte$X97.5.), ")", sep="")
tte$modelt <- paste(tte$model, "_t", tte$t, sep="")

tte$est <- tte$D
tte$D2.5 <- tte$X2.5.
tte$D97.5 <- tte$X97.5.
tte$grid <- tte$g
tte$model <- tte$modelt

tte.wide <- pivot_wider(tte, id_cols=c('sp','g'), names_from='modelt', values_from='d.ci')
tte.wide[tte.wide$sp == "TOWNSENDS CHIPMUNK",]$sp <- "chipmunk"
tte.wide[tte.wide$sp == "FLYING SQUIRREL",]$sp <- "flying\nsquirrel"
tte.wide[tte.wide$sp == "DEER MOUSE",]$sp <- "mouse"
tte.wide$sp <- tolower(tte.wide$sp)

###########################
sc.cam <- read.table(file="Output_cam/output_CAM-SC_separate_estimates.txt", sep=",", header=T) #separate
sc.cam2 <- read.table(file="Output_cam/output_CAM-SC_multigrid_estimates.txt", sep=",", header=T) #multigrid

#match sc.cam2 columns to sc.cam
sc.cam2$model <- sc.cam2$model.y
sc.cam2$mode <- sc.cam2$est

sc.cam <- rbind(sc.cam, sc.cam2[,names(sc.cam2) %in% names(sc.cam)])
sc.cam <- sc.cam[grep(sc.cam$param, pattern="D"),]

sc.cam <- sc.cam[!duplicated(sc.cam),] #check for and remove duplicates
sc.cam <- sc.cam[!is.na(sc.cam$Rhat),] #remove rows with NA for rhat
sc.cam <- sc.cam[sc.cam$Rhat < 1.2,] #remove rows with high rhat values
sc.cam$t <- factor(sc.cam$t, levels=c(0,15,60,1440))

# sc.cam$modelt <- paste("sc_cam ", sc.cam$model, "_t", sc.cam$t, sep="")
sc.cam$modelt <- paste(sc.cam$model, "_t", sc.cam$t, sep="")
sc.cam$mode.ci <- paste(round(sc.cam$mode, digits=2), " (", sprintf('%.2f', sc.cam$X2.5.), " - ", sprintf('%.2f', sc.cam$X97.5.), ")", sep="")

sc.cam$est <- sc.cam$mode
sc.cam$D2.5 <- sc.cam$X2.5.
sc.cam$D97.5 <- sc.cam$X97.5.
sc.cam$grid <- sc.cam$g
sc.cam$model <- sc.cam$modelt
# sc.cam[sc.cam$model == "sc_cam CAM-SC_info_multigrid_t0",]$model <- "SC-CAM_info_pool_t0"
sc.cam[sc.cam$model == "CAM-SC_info_multigrid_t0",]$model <- "SC-CAM_info_pool_t0"

sc.cam[grep(sc.cam$model, pattern="info_"),]$model <- gsub(sc.cam[grep(sc.cam$model, pattern="info_"),]$model, pattern="info_", replacement="info-")

sc.cam.wide <- pivot_wider(sc.cam, id_cols=c('sp','g'), names_from='modelt',values_from='mode.ci')


############################
#combine all into one data frame

unmarked <- merge(avgdets.wide, nmix.wide, by=c('sp','grid'), all=T)
unmarked <- merge(sp.convert, unmarked, by.x="acro", by.y="sp", all.y=T)
unmarked <- merge(unmarked, tte.wide, by.x=c('sp','grid'), by.y=c('sp','g'), all=T)
unmarked <- merge(unmarked, sc.cam.wide, by.x=c('acro','grid'), by.y=c('sp','g'), all=T)

unmarked[unmarked$sp == "flying\nsquirrel",]$sp <- "flying squirrel"

# write.table(unmarked, file="Data_linearregs/estimates_allunmarkedmodels.txt", sep=",", row.names = F)
############################

# u <- unmarked[,c("sp","grid","t0_d1","nmix pdecay_t15","tte_t1440","ste_t60","sc_cam CAM-SC_info_multigrid_t0")]
# u <- separate(u, col="t0_d1", into=c('avg.dets','sd'), sep=" \u00B1 ")
# u <- separate(u, col='nmix pdecay_t15', into=c('nmix pdecay_t15','nmix lwr','-','nmix upr'), sep=" ")
# u <- separate(u, col='tte_t1440', into=c('tte_t1440','nmix lwr','-','nmix upr'), sep=" ")

cols <- c('sp','grid','est','D2.5','D97.5','model')

u <- rbind(avgdets[avgdets$t_d == "avgdets_t0_d1", cols],
           nmix[nmix$mod_t == "nmix_pdecay_t15", cols],
           tte[tte$modelt %in% c("tte_t1440","ste_t60"),cols],
           sc.cam[sc.cam$modelt == "CAM-SC_info_multigrid_t0",cols])

u[u$sp == "PEMA",]$sp <- "DEER MOUSE"
u[u$sp == "TATO",]$sp <- "CHIPMUNK"
u[u$sp == "TOWNSENDS CHIPMUNK",]$sp <- "CHIPMUNK"
u[u$sp == "GLOR",]$sp <- "FLYING SQUIRREL"
u[u$sp == "GLSA",]$sp <- "FLYING SQUIRREL"

random <- read.table(file="Data_linearregs/CR_estimates_correctgridnums.txt", sep=",", header=T)
random <- random[random$model == "SCR random",] #has absolute values and scaled values
random <- random[,c("sp","grid",
                    "dhat","mean.random","sd.random","X2.5..random","X97.5..random", #remove extra columns
                    "dhat.scaled","X2.5.r.scaled","X97.5.r.scaled")]
random <- random[,c('sp','grid','dhat','X2.5..random','X97.5..random')]
random$model <- "SCR_random"
random$sp <- toupper(random$sp)
random[random$sp == "MOUSE",]$sp <- "DEER MOUSE"
random[random$sp == "FLYINGSQUIRREL",]$sp <- "FLYING SQUIRREL"
names(random) <- cols

u <- rbind(u, random)
u$model <- factor(u$model, 
                  levels=c("avgdets_t0_d1","nmix_pdecay_t15","tte_t1440","ste_t60","SC-CAM_info-pool_t0","SCR_random"),
                  labels=c("avgdets_t1_d1","nmix_pdecay_t15","tte_t1440","ste_t60","SC-CAM_info-pool_t1","SCR_random"))
u$sp <- factor(u$sp, levels=c("DEER MOUSE","CHIPMUNK","FLYING SQUIRREL"))

ggplot(u[!u$model %in% "SCR_random",], aes(x=grid, y=est, col=model, group=model)) + 
  geom_point(size=3, position=position_dodge(width=0.5)) +   geom_errorbar(aes(ymin=D2.5, ymax=D97.5), width=0.1, position=position_dodge(width=0.5)) + 
  geom_point(data=u[u$model == "SCR_random",], col="black", size=3) + geom_errorbar(data=u[u$model == "SCR_random",], aes(ymin=D2.5, ymax=D97.5), width=0.1, position=position_dodge(width=0.5), col="black") + 
  scale_x_continuous(breaks=1:9, name="site") +
  scale_y_continuous(name="model estimate") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="top") + facet_wrap(~sp, scales = "free_x") + coord_flip()
# ggsave("Figures/unmarked_estimates.tiff", height=6, width=12, units="in", compression="lzw", dpi=400)

#######################################
#calculate Coefficient of variation (CV)

#multiple definitions of CV:
#1. relative standard deviation (RSD)
#2. relative standard error (RSE)

#Proportional Standard Error as a measure of precision
# CV(theta_hat) = SE(theta_hat)/theta, standard error of the estimator of theta divided by theta, the parameter itself
#Pollock et al. 1990
#lowering CV indicates increasing precision for the study. Rough rule of thumb is that a study that provides a CV of 20% is reasonable.

avgdets$CV <- avgdets$se/avgdets$mean * 100
avgdets$metric <- "RSE"

nmix$CV <- nmix$sd/nmix$mean * 100
nmix$metric <- "RSD"

tte$CV <- tte$SE/tte$N * 100
tte$metric <- "RSE"

sc.cam$CV <- sc.cam$sd/sc.cam$mean * 100
sc.cam$metric <- "RSD"

cols <- c("sp","grid","model","CV","metric")

cv <- rbind(avgdets[,cols], nmix[,cols], tte[,cols], sc.cam[,cols])
  cv[cv$sp %in% "PEMA",]$sp <- "DEER MOUSE"
  cv[cv$sp %in% "TATO",]$sp <- "TOWNSENDS CHIPMUNK"
  cv[cv$sp %in% "GLSA",]$sp <- "FLYING SQUIRREL"
  cv[cv$sp %in% "GLOR",]$sp <- "FLYING SQUIRREL"
  
  cv <- separate(cv, col="model", remove=F, sep="_", into=c("mod","type","t"))
  
  cv[cv$mod == "tte",]$t <- cv[cv$mod == "tte",]$type
  cv[cv$mod == "tte",]$type <- ""
  
  cv[cv$mod == "ste",]$t <- cv[cv$mod == "ste",]$type
  cv[cv$mod == "ste",]$type <- ""
  
  cv$d <- ""
  cv[cv$mod == "avgdets",]$d <- cv[cv$mod == "avgdets",]$t
  cv[cv$mod == "avgdets",]$t <- cv[cv$mod == "avgdets",]$type
  cv[cv$mod == "avgdets",]$type <- cv[cv$mod == "avgdets",]$d
  
  cv[cv$mod == "CAM-SC",]$mod <- "SC-CAM"
  
  cv[cv$t == "t0",]$t <- "t1"
  
  cv[grep(cv$type, pattern="multigrid"),]$type <- gsub(cv[grep(cv$type, pattern="multigrid"),]$type, pattern="multigrid", replacement="pool")
  
  #add NA values
  # cv.empty1 <- data.frame(sp=rep(c("DEER MOUSE","TOWNSENDS CHIPMUNK","FLYING SQUIRREL"), times=3), 
  #                        t=c(rep("t1440", times=3), rep("t60", times=3), rep("t15", times=3)), # , rep("t1", times=3)), 
  #                        type=rep("info\npool", times=9), 
  #                        grid=0, CV=-10, d="d0", mod="SC-CAM", model="empty", metric=NA)
  cv.empty1 <- data.frame(sp=rep(c("FLYING SQUIRREL"), times=3), 
                          t=c("t15","t60","t1440"), 
                          type=rep("info\npool", times=3), 
                          grid=0, CV=-10, d="d0", mod="SC-CAM", model="empty", metric=NA)
  cv.empty2 <- data.frame(sp=rep(c("DEER MOUSE","TOWNSENDS CHIPMUNK","FLYING SQUIRREL"), times=3), 
                         t=c(rep("t1", times=3), rep("t60", times=3), rep("t15", times=3)), # , rep("t1", times=3)), 
                         type=rep("info\npool", times=9), 
                         grid=0, CV=-10, d="d0", mod="SC-CAM", model="empty", metric=NA)
  cv.empty2 <- cv.empty2[!(cv.empty2$sp == "DEER MOUSE" & cv.empty2$t %in% c("t15","t60","t1440")),] #remove those with actual estimates
  cv.empty3 <- data.frame(sp="FLYING SQUIRREL", 
                          t="t1", 
                          type="noinfo\npool", 
                          grid=0, CV=-10, d="d0", mod="SC-CAM", model="empty", metric=NA)
  cv.empty4 <- data.frame(sp="TOWNSENDS CHIPMUNK", 
                          t="t1440", 
                          type="noinfo\nseparate", 
                          grid=0, CV=-10, d="d0", mod="SC-CAM", model="empty", metric=NA)
  cv.empty5 <- data.frame(sp="DEER MOUSE", 
                          t="t1", 
                          type="info\nseparate", 
                          grid=0, CV=-10, d="d0", mod="SC-CAM", model="empty", metric=NA)
  
  cv <- rbind(cv, cv.empty1, cv.empty2, cv.empty3, cv.empty4, cv.empty5)
  
  table(cv[,c("mod","t")])
  
  #factor variables
  cv$mod <- factor(cv$mod, levels=c("avgdets","nmix","tte","ste","SC-CAM"))
  cv$t <- factor(cv$t, levels=c("t1","t5","t15","t30","t60","t1440"))
  cv$sp <- factor(cv$sp, levels=c("DEER MOUSE","TOWNSENDS CHIPMUNK","FLYING SQUIRREL"))
  
  cv$type <- gsub(cv$type, pattern="-", replacement="\n")
  cv$type <- factor(cv$type, levels=c(paste("d",1:9, sep=""), "base","pstation","pdecay","","noinfo\nseparate","noinfo\npool","info\nseparate","info\npool"))

# ddply(cv, .(sp, model, metric), summarize, cv.mean=mean(CV, na.rm=T))

ggplot(data=cv) + 
  geom_hline(aes(yintercept = 0), lwd=0.1, col="black") + 
  geom_hline(aes(yintercept = 100), lwd=0.1, col="black") + 
  geom_hline(aes(yintercept = 20), lty="dashed", lwd=0.5, col="grey50") + 
  geom_boxplot(aes(x=type, y=CV, fill=t), width=0.75, outliers = F) + 
  geom_point(aes(x=type, y=CV, fill=t, group=t), position=position_jitterdodge(jitter.width = 0.3), pch=21, col="black", alpha=0.5) +
  
  # geom_boxplot(data=cv.empty, aes(x=type, y=CV, group=t), width=0.75, outliers = F, fill=NA) + 
  facet_grid(cols=vars(mod), rows=vars(sp), scales="free_x", space="free_x") +
  # coord_cartesian(ylim=c(0,100)) +
  coord_cartesian(ylim=c(0, NA)) +
  # facet_wrap(~sp + mod, scales="free_x", nrow=3) + 
  ylab("Precision (%)") + xlab("") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank())
ggsave("Figures/CV_unmarked.tiff", height=12, width=25, units="in", compression="lzw", dpi=350)
