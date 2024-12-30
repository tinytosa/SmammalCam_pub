
#########################
#combine all data from capture-recapture models
#need to plot on log-log scale

#remember that no shermans and PEMA on grid 6. Need to fix grid numbers for grid 7, 8, and 9

#########################
#load packages
require(xlsx)
require(dplyr) #for row_binds
require(tidyr)

require(ggplot2)
require(ggpubr)

#########################
c.vals <- c("black","grey70", #"steelblue","#83B6C3",
                   "#6F1A07","#EE2677","#FCBFB7",
                   "#3D523D","#5F815F","#AAC0AA")
                   

# order.mods <- c("MNKA","Huggins","SCR_separate","SCR_allgrids","SCR_random","SC noinfo","SC info","SC noinfo all","SC info all")
# order.mods.names <- c("MNKA","Huggins","SCR separate","SCR pool","SCR random","SC no.info sep","SC info sep","SC no.info pool","SC info pool")

#remove SC models
order.mods <- c("MNKA","Huggins","SCR_separate","SCR_allgrids","SCR_random")
order.mods.names <- c("MNKA","Huggins","SCR separate","SCR pool","SCR random")

order.sp <- c("PEMA","TATO","GLSA","GLSA")
order.sp.names <- c("mouse","chipmunk","flyingsquirrel","flying\nsquirrel")

#########################
#skip to below
#########################
#load data

mnka <- read.table(file="C:/Users/tosam/Documents/0_OSU/Dissertation/Data/SmallMammalCams_SCR/Output_SCR/MNKA.txt", sep=",") #PEMA grid numbers are incorrect

huggins <- read.xlsx("C:/Users/tosam/Documents/0_OSU/Dissertation/Data/SmallMammalCams/SmallMammal_Huggins_Estimates_fromWeldy_updated.xlsx", sheetName="Abundances", colClasses="character") #PEMA grid numbers are correct
#huggins is missing grid 6 for PEMA

separate <- read.table(file="C:/Users/tosam/Documents/0_OSU/Dissertation/Data/SmallMammalCams_SCR/Output_separate/SCR-B_separate.txt", sep=",") #grid numbers for PEMA are correct

allgrids.mean <- read.table(file="C:/Users/tosam/Documents/0_OSU/Dissertation/Data/SmallMammalCams_SCR/Output_SCR/SCR-B_allgrids.txt", sep=",")
allgrids.mode <- read.table(file="C:/Users/tosam/Documents/0_OSU/Dissertation/Data/SmallMammalCams_SCR/Output_SCR/SCR-B_allgrids_mode.txt", sep=",")

random.mean <- read.table(file="C:/Users/tosam/Documents/0_OSU/Dissertation/Data/SmallMammalCams_SCR/Output_SCR/SCR-B_random.txt", sep=",")
random.mode <- read.table(file="C:/Users/tosam/Documents/0_OSU/Dissertation/Data/SmallMammalCams_SCR/Output_SCR/SCR-B_random_mode.txt", sep=",")

#########################
#calculate CV?
# Coefficient of variation (CV) = relative standard deviation (RSD)
# a measure of precision where lower values indicates greater precision
# CV = sd/mu or se/mu
# CI = mu +/- z*s/sqrt(n)

#Relative standard error (RSD)
#RSE = SE(D)/D
#RSE can be equivalent to %CV

#can't for mnka

#huggins
huggins$CV <- huggins$SE/huggins$N * 100
huggins$metric <- "RSE"
huggins$sp <- "DEER MOUSE"
huggins[huggins$Sp == "TATO",]$sp <- "TOWNSENDS CHIPMUNK"
huggins[huggins$Sp == "GLSA",]$sp <- "FLYING SQUIRREL"
huggins$grid <- huggins$Grid
huggins$model <- "Huggins"
huggins <- huggins[grep(huggins$ModelName, pattern="Null"),]

separate$CV <- separate$sd/separate$mean * 100
separate$metric <- "RSD"
separate <- separate[separate$param == "D",]
separate[separate$sp == "PEMA",]$sp <- "DEER MOUSE"
separate[separate$sp == "TATO",]$sp <- "TOWNSENDS CHIPMUNK"
separate[separate$sp == "GLSA",]$sp <- "FLYING SQUIRREL"
separate$model <- "SCR_separate"

allgrids.mean$CV <- allgrids.mean$sd/allgrids.mean$mean * 100
allgrids.mean$metric <- "RSD"
allgrids.mean <- allgrids.mean[grep(allgrids.mean$param, pattern="D"),]
allgrids.mean[allgrids.mean$sp == "PEMA",]$sp <- "DEER MOUSE"
allgrids.mean[allgrids.mean$sp == "TATO",]$sp <- "TOWNSENDS CHIPMUNK"
allgrids.mean[allgrids.mean$sp == "GLSA",]$sp <- "FLYING SQUIRREL"
allgrids.mean$model <- "SCR_pool"

random.mean$CV <- random.mean$sd/random.mean$mean * 100
random.mean$metric <- "RSD"
random.mean <- random.mean[grep(random.mean$param, pattern="D"),]
random.mean[random.mean$sp == "PEMA",]$sp <- "DEER MOUSE"
random.mean[random.mean$sp == "TATO",]$sp <- "TOWNSENDS CHIPMUNK"
random.mean[random.mean$sp == "GLSA",]$sp <- "FLYING SQUIRREL"

cols <- c("sp","model","CV","metric")

cv <- rbind(huggins[,cols], separate[,cols], allgrids.mean[,cols], random.mean[,cols])
cv$model <- gsub(cv$model, pattern="_", replacement="\n")
cv$model <- factor(cv$model, levels=c("Huggins","SCR\nseparate","SCR\npool","SCR\nrandom"))
cv$sp <- factor(cv$sp, levels=c("DEER MOUSE","TOWNSENDS CHIPMUNK","FLYING SQUIRREL"))

ggplot(data=cv, aes(x=model, y=CV, fill=model)) + 
  geom_hline(aes(yintercept = 0), lwd=0.1, col="black") + 
  # geom_hline(aes(yintercept = 100), lwd=0.1, col="black") + 
  geom_hline(aes(yintercept = 20), lty="dashed", lwd=0.5, col="grey50") + 
  
  geom_boxplot(outliers = F) + 
  geom_point(aes(fill=model), position=position_jitterdodge(jitter.width = 0.3),  pch=21, col="black", alpha=0.5) + 
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position = "", strip.text.x = element_blank()) +
  # coord_cartesian(ylim=c(0,100)) +
  # coord_cartesian(ylim=c(0,75)) +
  ylab("Precision (%)") + xlab("") +
  facet_grid(col=vars(model), row=vars(sp), scales="free_x", space="free_x")
# ggsave("Figures/CV_mark-recapture.tiff", height=10, width=6, units="in", compression="lzw", dpi=350)

cv[cv$model == "SCR\nrandom",]

#########################
#calculate area of inference using mmdm and 1/2 mmdm

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

#########################
#columns to keep
c2k <- c("sp","grid","model","dhat.random","mean.random","sd.random","X2.5..random","X50..random","X97.5..random","dhat","X2.5.","X97.5.","Rhat")

random.mean$grid <- substr(random.mean$param, start=nchar(random.mean$param)-1, stop=nchar(random.mean$param)-1) #extract grid info
random.mean <- random.mean[grep(random.mean$param, pattern="D"),] #subset for density
random.mean <- random.mean[,c("mean","sd","X2.5.","X50.","X97.5.","Rhat","sp","grid","model")]
names(random.mean) <- c("mean.random","sd.random","X2.5..random","X50..random","X97.5..random","Rhat.random","sp","grid","model.random")

random.mode <- random.mode[,c("dhat","nhat","sp","grid")]
names(random.mode) <- c("dhat.random","nhat.random","sp","grid")

random <- merge(random.mode, random.mean, by=c("sp","grid"), all.x=T)

random2 <- merge(random, random, by=c("sp","grid"), all.x=T, suffixes=c(".og",""))
names(random2) <- gsub(names(random2), pattern=".random.og", replacement="")

#
random %>% dplyr::group_by(sp) %>%
  reframe(dhat.range=range(dhat.random), nhat.range=range(nhat.random)) %>% as.data.frame() %>% mutate_if(is.numeric, round, digits=2)

#     sp dhat.range nhat.range
# 1 GLSA       1.14      63.94
# 2 GLSA       3.02     168.82
# 3 PEMA      13.12      65.92
# 4 PEMA      62.84     315.64
# 5 TATO       1.68      88.00
# 6 TATO       8.59     448.72
#
mnka <- merge(mnka, random, by.x=c("sp","g"), by.y=c("sp","grid"), all.x=T)
mnka <- merge(mnka, info[, c("sp","mmdm.a","halfmmdm.a","n.occ.ideal")], by="sp", all.x=T)
mnka$dhat <- mnka$mnka/mnka$mmdm.a #convert from n to density
# mnka$dhat <- mnka$mnka/mnka$mmdm.a/mnka$noccasions*mnka$n.occ.ideal#account for fewer occassions for some grids
# mnka$dhat <- mnka$mnka/(280*280/10000) #convert N to density, area = 280*280 m^2, 10000 m^2 = 1 ha for large grid
# mnka[mnka$sp == "PEMA",]$mnka.d <- mnka[mnka$sp == "PEMA",]$mnka/(90*90/10000)
names(mnka) <- gsub(names(mnka), pattern="g", replacement="grid")
mnka$model <- "MNKA"

#
huggins <- huggins[huggins$ModelName == "Huggins Null",] #filter out huggins top models
names(huggins) <- c("Species","sp","mod","model","grid","n","sd","n.ci2.5","n.ci97.5","notes") #change
# huggins$g <- huggins$grid
# huggins$grid <- paste("SM-0", huggins$grid, sep="")
huggins$sp <- as.character(huggins$sp)
# huggins$mod <- as.character(huggins$mod)
huggins$model <- "Huggins"

huggins <- merge(huggins, info[,c("sp","mmdm.a","halfmmdm.a")], by="sp", all.x=T)
huggins$dhat <- huggins$n/huggins$mmdm.a
huggins$X2.5. <- huggins$n.ci2.5/huggins$mmdm.a
huggins$X97.5. <- huggins$n.ci97.5/huggins$mmdm.a

#huggins grid numbers are correct. 
huggins[huggins$sp == "PEMA",]$grid <- 1:7
huggins <- merge(huggins, random, by=c("sp","grid"), all.x=T)

#
#grid numbers for PEMA are correct
separate <- separate[separate$param == "D",]
separate$model <- "SCR_separate"
#for merging with random, need to make grid numbers for PEMA match the ones in random, easiest thing to do is turn separate grid numbers to incorrect ones
separate[separate$sp == "PEMA",]$grid <- 1:8
separate <- merge(separate, random, by=c("sp","grid"), all.x=T)
separate$dhat <- separate$mode

#
allgrids.mean$grid <- substr(allgrids.mean$param, start=nchar(allgrids.mean$param)-1, stop=nchar(allgrids.mean$param)-1)
allgrids.mean <- allgrids.mean[grep(allgrids.mean$param, pattern="D"),] #subset for density
allgrids.mean <- allgrids.mean[,c("mean","sd","X2.5.","X50.","X97.5.","Rhat","sp","grid","model")]

allgrids.mode <- allgrids.mode[,c("dhat","nhat","sp","grid")]
# names(allgrids.mode) <- c("dhat.allgrids","nhat.allgrids","sp","grid")

allgrids <- merge(allgrids.mode, allgrids.mean, by=c("sp","grid"), all.x=T)
allgrids <- merge(allgrids, random, by=c("sp","grid"), all.x=T)

############
#combine into one data frame
allscr <- bind_rows(random2[,names(random2) %in% c2k], allgrids[,names(allgrids) %in% c2k], separate[,names(separate) %in% c2k], huggins[,names(huggins) %in% c2k], mnka[,names(mnka) %in% c2k]) #, unmcr[,names(unmcr) %in% c2k])
table(allscr$model)

allscr$converged <- "converged"
try(allscr[allscr$Rhat > 1.1 & !is.na(allscr$Rhat),]$converged <- "not converged")

allscr <- allscr[allscr$converged == "converged",] #remove non-converged 

############
#create table for SI with all CR values

#create range column
allscr$range <- paste(round(allscr$X2.5., digits=2), "-", round(allscr$X97.5., digits=2), sep=" ")
allscr[allscr$model == "MNKA",]$range <- ""

allscr$est_range <- paste(round(allscr$dhat, digits=2), " (", allscr$range, ")", sep="")
allscr[allscr$model == "MNKA",]$est_range <- round(allscr[allscr$model == "MNKA",]$dhat, digits=2)

#convert from long to wide
allscr.wide <- pivot_wider(data=allscr[,c("sp","grid","model","est_range")], names_from=c(model), values_from=est_range)
allscr.wide[allscr.wide$sp == "PEMA",]$grid <- c(1:5,7:9) #correct wrong grid numbers for PEMA

allscr.wide <- allscr.wide[,c("sp","grid","MNKA","Huggins","SCR_separate","SCR_allgrids","SCR_random")]

# write.csv(allscr.wide, file="Output_CR/all_CR.csv")

############
#1. plot
random$sp <- factor(random$sp, levels=c("PEMA","TATO","GLSA","GLOR"))
random[random$sp == "PEMA",]$grid <- c(1:5,7:9) #correct wrong grid numbers for PEMA
ggplot(data=random, aes(x=grid, y=dhat.random)) + geom_errorbar(aes(ymin=X2.5..random, ymax=X97.5..random), width=0.5) + geom_point(size=3) + theme_bw(base_size = 20) + facet_wrap(~sp, scales="free_y") + scale_y_log10() + scale_x_continuous(breaks=c(1:9)) + theme(panel.grid = element_blank()) #all SCR random estimates by grid

#
# sp.cr <- "PEMA"
# ggplot() +
#   geom_errorbar(data=random2[random2$sp == sp.cr,], aes(x=X50..random, ymin=X2.5..random, ymax=X97.5..random), col="#6B6D76", linewidth=1) + 
#   geom_errorbar(data=random2[random2$sp == sp.cr,], aes(y=dhat.random, xmin=X2.5., xmax=X97.5.), col="#6B6D76", linewidth=1) + 
#   geom_point(data=random2[random2$sp == sp.cr,], aes(x=X50..random, y=dhat.random), size=3, col="#6B6D76") + 
#   geom_smooth(data=random2[random2$sp == sp.cr,], aes(x=X50..random, y=dhat.random), method="lm", alpha=0.1) + #SCRB-mode vs SCRB-median for random effect
#   # geom_errorbar(data=allgrids[allgrids$sp == sp.cr,], aes(x=dhat, ymin=X2.5..random, ymax=X97.5..random), col="#334E58", linewidth=1, width=1) + geom_point(data=allgrids[allgrids$sp == sp.cr,], aes(x=dhat, y=dhat.random), size=3, col="#334E58") + #SCRB-allgrids vs SCRBrandom
#   # geom_errorbar(data=mnka[mnka$sp == sp.cr,], aes(x=dhat, ymin=X2.5..random, ymax=X97.5..random), col="#410200", linewidth=1) + geom_point(data=mnka[mnka$sp == sp.cr,], aes(x=dhat, y=dhat.random), size=3, col="#410200") + #mnka vs SCRBrandom
#   #huggins vs SCRBrandom
#   geom_abline(slope=1, intercept=0, col="grey60", lty="dashed") +
#   xlab("model estimate") + ylab("SCR random (animals/ha)") +
#   theme_bw(base_size=20) + theme(panel.grid = element_blank())

allscr$model <- factor(allscr$model, levels=order.mods, labels=order.mods.names)
allscr$sp <- factor(allscr$sp, levels=order.sp, labels=order.sp.names)

plot.cr <- function(sp.cr)
{
  plot <- ggplot(data=allscr[allscr$sp == sp.cr & !allscr$model %in% c("SC no.info","SC info","SC info pool"),], aes(x=dhat, y=dhat.random, col=model)) +
    geom_abline(slope=1, intercept=0, col="grey60", linewidth=0.25) +
    geom_errorbar(aes(ymin=X2.5..random, ymax=X97.5..random), width=0.01, lty="dashed", linewidth=0.1) + #vertical error bars
    geom_errorbar(aes(xmin=X2.5., xmax=X97.5.), width=0.01, lty="dashed", linewidth=0.1) + #horizontal error bars
    geom_point(size=3) + #, alpha=0.5) +
    geom_smooth(method="lm", aes(fill=model, lty=model), alpha=0.1) + #add ribbons
    scale_color_manual(values=c.vals) + scale_fill_manual(values=c.vals) +
    # scale_x_log10() + scale_y_log10() +
    # scale_x_sqrt() + scale_y_sqrt() +
    xlab("model estimate (animals/ha)") + ylab("SCR random (animals/ha)") +
    theme_bw(base_size=20) + theme(panel.grid = element_blank()) + coord_equal()
  # ggsave(plot, filename = paste("Figures/CR_density_", sp.cr, ".tiff"), height=5, width=8, units="in", dpi=400, compression="lzw")
  return(plot)
}
# p <- plot.cr("PEMA")
# t <- plot.cr("TATO")
# g <- plot.cr("GLSA")
p <- plot.cr("mouse")
t <- plot.cr("flying\nsquirrel")
g <- plot.cr("chipmunk")

#code below works better
# all <- ggarrange(p, t+rremove("ylab"), g+rremove("ylab"), nrow=1, align="hv", common.legend = T, legend = "top") #, labels="AUTO", font.label = list(size=32))
# all
# ggsave(all, filename="Figures/CR_density.tiff", height=6, width=20, units="in", dpi=400, compression="lzw")

all <-  ggplot(data=allscr[!allscr$model %in% c("SC no.info","SC info","SC info pool"),] , aes(x=dhat, y=dhat.random, col=model)) +
  geom_abline(slope=1, intercept=0, col="grey60", linewidth=0.25) +
  geom_errorbar(aes(ymin=X2.5..random, ymax=X97.5..random), width=0.01, lty="dashed", linewidth=0.1) + #vertical error bars
  geom_errorbar(aes(xmin=X2.5., xmax=X97.5.), width=0.01, lty="dashed", linewidth=0.1) + #horizontal error bars
  geom_point(size=3, pch=1, stroke=2) + #, alpha=0.5) +
  geom_smooth(method="lm", aes(fill=model, lty=model), alpha=0.1) + #add ribbons
  scale_color_manual(values=c.vals) + scale_fill_manual(values=c.vals) +
  # scale_x_log10() + scale_y_log10() +
  # scale_x_sqrt() + scale_y_sqrt() +
  xlab("model estimate (animals/ha)") + ylab("SCR random (animals/ha)") +
  theme_bw(base_size=20) + theme(panel.grid = element_blank(), legend.position="top", strip.background = element_blank(), strip.text.x = element_blank()) + facet_wrap(~sp, scales="free")
all
# ggsave(all, filename="Figures/CR_density.tiff", height=6, width=15, units="in", dpi=400, compression="lzw")

#for CR-SC models, but do this later in another script. ignore for now
# SCall <-  ggplot(data=allscr[allscr$model %in% c("SC no.info","SC info","SC info pool"),] , aes(x=dhat, y=dhat.random, col=model)) +
#   geom_abline(slope=1, intercept=0, col="grey60", linewidth=0.25) +
#   geom_errorbar(aes(ymin=X2.5..random, ymax=X97.5..random), width=0.01, lty="dashed", linewidth=0.1) + #vertical error bars
#   geom_errorbar(aes(xmin=X2.5., xmax=X97.5.), width=0.01, lty="dashed", linewidth=0.1) + #horizontal error bars
#   geom_point(size=3, pch=1, stroke=2) + #, alpha=0.5) +
#   geom_smooth(method="lm", aes(fill=model, lty=model), alpha=0.1) + #add ribbons
#   scale_color_manual(values=sc.vals) + scale_fill_manual(values=sc.vals) +
#   # scale_x_log10() + scale_y_log10() +
#   # scale_x_sqrt() + scale_y_sqrt() +
#   xlab("model estimate (animals/ha)") + ylab("SCR random (animals/ha)") +
#   theme_bw(base_size=20) + theme(panel.grid = element_blank(), legend.position="top", strip.background = element_blank(), strip.text.x = element_blank()) + facet_wrap(~sp, scales="free")
# SCall
# ggsave(SCall, filename="Figures/CR_SC_density.tiff", height=6, width=15, units="in", dpi=400, compression="lzw")

############################
#2. calculate slopes

#linear regression assumptions:
#1. linearity of data
#2. normality of residuals
#3. homogeneity of resid variance
#4. independence of residual error terms

#square root transform or log-log transform and fit linear regression

unique(allscr$model)
# SCR random   SCR allgrids      SCR separate Huggins      MNKA

mod.coef.all <- NULL
mod.coef.scaled <- NULL
allscr.scaled <- NULL

# sp.cr <- "PEMA"
# sp.cr <- "GLSA"
# sp.cr <- "TATO"
# mod <- "SCR all"

###run linear regressions for each species for each model
# for(sp.cr in c("PEMA","GLSA","TATO"))
for(sp.cr in c("mouse","flying\nsquirrel","flyingsquirrel","chipmunk"))
{
  for(mod in unique(allscr$model))
  # for(mod in unique(allscr$model)[2:5]) #after TATO SCR random
  {
    print(paste(sp.cr, mod))
    # if(sp.cr == "flyingsquirrel") #since errors with returns in sp.cr name
    # { 
    #   sp.cr <- "flying\nsquirrel"
    # }
    
    data <- allscr[allscr$sp == sp.cr & allscr$model == mod,]
    if(nrow(data) == 0)
    {
      next
    }
    if(sp.cr == "flying\nsquirrel") #since errors with returns in sp.cr name
    { 
      sp.cr <- "flyingsquirrel"
    }
    l <- lm(dhat.random ~ dhat, data=data)
    # l.sqrt <- lm(sqrt(dhat.random) ~ sqrt(dhat), data=data) #don't need to do since the residuals don't look too bad
    # l.log <- lm(log(dhat.random) ~ log(dhat), data=data)
    
    #scaled variables
    data$dhat.scaled <- scale(data$dhat) #center and scale independent var
    data$X2.5.scaled <- (data$X2.5. - attr(data$dhat.scaled, "scaled:center"))/attr(data$dhat.scaled, "scaled:scale")
    data$X97.5.scaled <- (data$X97.5. - attr(data$dhat.scaled, "scaled:center"))/attr(data$dhat.scaled, "scaled:scale")
    
    data$dhat.random.scaled <- scale(data$dhat.random) #center and scale dependent var
    data$X2.5.r.scaled <- (data$X2.5..random - attr(data$dhat.random.scaled, "scaled:center"))/attr(data$dhat.random.scaled, "scaled:scale")
    data$X97.5.r.scaled <- (data$X97.5..random - attr(data$dhat.random.scaled, "scaled:center"))/attr(data$dhat.random.scaled, "scaled:scale")
    
    allscr.scaled <- rbind(allscr.scaled, data)
    
    if(nrow(data) == 1)
    {
      next
    }
      
    if(nrow(data) > 2)
    {
      #look at residuals using plot
      tiff(filename=paste("Figures/residuals_plots/", mod, "_", sp.cr,".tiff", sep=""), height=8, width=8, units="in", res=300, compression="lzw")
      par(mfrow=c(2,2))
      plot(l)
      # plot(l.sqrt)
      # plot(l.log)
      dev.off()
      
      #look at residuals using ggplot
      # require(ggfortify) 
      # autoplot(l) + theme_bw(base_size=20) + theme(panel.grid = element_blank())
      # p <- autoplot(l.sqrt) + theme_bw(base_size=20) + theme(panel.grid = element_blank())
      # autoplot(l.log) + theme_bw(base_size=20) + theme(panel.grid = element_blank())
      # ggsave(p, filename=paste("Figures/residuals_plots/", mod, "_", sp.cr,".tiff", sep=""), height=8, width=8, units="in", dpi=300, compression="lzw")
    }
    
    mod.coef <- data.frame(coef(summary(l)))[2,]
    mod.coef$sp <- sp.cr
    mod.coef$mod <- mod
    mod.coef$r2 <- summary(l)$adj.r.squared
    mod.coef$rmse <- sqrt(mean(l$residuals^2))
    mod.coef.all <- rbind(mod.coef.all, mod.coef)
    
    #scaled linear regression
    l <- lm(dhat.random.scaled ~ dhat.scaled, data=data) #linear regression
    if(nrow(data) > 2)
    {
      #look at residuals using plot
      tiff(filename=paste("Figures/residuals_plots/", mod, "_", sp.cr,"_scaled.tiff", sep=""), height=8, width=8, units="in", res=300, compression="lzw")
      par(mfrow=c(2,2))
      plot(l)
      dev.off()
    }
    
    mod.coef.s <- data.frame(coef(summary(l)))[2,]
    mod.coef.s$sp <- sp.cr
    mod.coef.s$mod <- mod
    mod.coef.s$r2 <- summary(l)$adj.r.squared
    mod.coef.s$rmse <- sqrt(mean(l$residuals^2))
    mod.coef.scaled <- rbind(mod.coef.scaled, mod.coef.s)
  }
}

#fix dhat.scaled and dhat.random.scaled when only 1 data point for sp.cr and mod
# allscr.scaled[is.nan(allscr.scaled$dhat.scaled),]$dhat.scaled <- 0
# allscr.scaled[is.nan(allscr.scaled$dhat.random.scaled),]$dhat.random.scaled <- 0

allscr.scaled <- allscr.scaled[!duplicated(allscr.scaled),] #129 rows

# write.table(allscr.scaled, file="Data_linearregs/CR_estimates.txt", sep=",", row.names=F)

#correct PEMA grid numbers
allscr.scaled[allscr.scaled$sp == "mouse" & allscr.scaled$grid == 8,]$grid <- 9
allscr.scaled[allscr.scaled$sp == "mouse" & allscr.scaled$grid == 7,]$grid <- 8
allscr.scaled[allscr.scaled$sp == "mouse" & allscr.scaled$grid == 6,]$grid <- 7
# write.table(allscr.scaled, file="Data_linearregs/CR_estimates_correctgridnums.txt", sep=",", row.names=F)

#fix flying squirrel sp
mod.coef.all[mod.coef.all$sp == "flyingsquirrel",]$sp <- "flying\nsquirrel"

#warning for SCR random is correct since it should be essentially a perfect fit
mod.coef.all %>% mutate_if(is.numeric, round, digits=2)

#        Estimate Std..Error      t.value Pr...t..               sp          mod   r2  rmse
# dhat       1.00       0.00 2.214744e+16     0.00            mouse   SCR random 1.00  0.00
# dhat1      1.30       0.24 5.350000e+00     0.00            mouse     SCR pool 0.80  6.69
# dhat2      1.03       0.08 1.276000e+01     0.00            mouse SCR separate 0.96  3.03
# dhat3      1.55       0.28 5.510000e+00     0.00            mouse      Huggins 0.83  6.19
# dhat4      1.68       0.59 2.850000e+00     0.03            mouse         MNKA 0.50 10.47

# dhat5      1.00       0.00 1.340617e+16     0.00 flying\nsquirrel   SCR random 1.00  0.00
# dhat6      1.01       0.09 1.160000e+01     0.00 flying\nsquirrel     SCR pool 0.94  0.15
# dhat7      0.92       0.07 1.263000e+01     0.00 flying\nsquirrel SCR separate 0.95  0.13
# dhat8      1.33       0.13 1.056000e+01     0.00 flying\nsquirrel      Huggins 0.93  0.16
# dhat9      1.38       0.16 8.510000e+00     0.00 flying\nsquirrel         MNKA 0.90  0.19

# dhat10     1.00       0.00          Inf     0.00         chipmunk   SCR random 1.00  0.00
# dhat11     1.01       0.07 1.515000e+01     0.00         chipmunk     SCR pool 0.97  0.42
# dhat12     1.00       0.02 4.754000e+01     0.00         chipmunk SCR separate 1.00  0.14
# dhat13     1.27       0.10 1.300000e+01     0.00         chipmunk      Huggins 0.95  0.49
# dhat14     1.48       0.11 1.311000e+01     0.00         chipmunk         MNKA 0.96  0.48

mod.coef.all$mod <- factor(mod.coef.all$mod, levels=order.mods.names, labels=order.mods.names)
# mod.coef.all$sp <- factor(mod.coef.all$sp, levels=order.sp, labels=order.sp.names)
mod.coef.all$sp <- factor(mod.coef.all$sp, levels=order.sp.names, labels=order.sp.names)

# mod.coef.all.sc <- mod.coef.all[mod.coef.all$mod %in% c("SC info pool","SC no.info","SC info"),]
mod.coef.all <- mod.coef.all[!mod.coef.all$mod %in% c("SC info pool","SC no.info","SC info"),]

# write.table(mod.coef.all, file="Data_linearregs/CR_slopes_individualidentity.txt", sep=",", row.names=F)
# write.table(mod.coef.scaled, file="Data_linearregs/CR_slopes_scaled.txt", sep=",", row.names=F)

#removed this
# write.table(mod.coef.all.sc, file="Data_linearregs/CR_slopes_spatialcount.txt", sep=",", row.names=F)


#####
#summarize
#####
mod.coef.all %>% dplyr::group_by(sp) %>%
  reframe(range.slope=range(Estimate), range.r2=range(r2), range.rmse=range(rmse)) %>% as.data.frame() %>% mutate_if(is.numeric, round, digits=2)

#                 sp range.slope range.r2 range.rmse
# 1            mouse        1.00     0.50       0.00
# 2            mouse        1.68     1.00      10.47

# 3         chipmunk        1.00     0.95       0.00
# 4         chipmunk        1.48     1.00       0.49

# 5 flying\nsquirrel        0.92     0.90       0.00
# 6 flying\nsquirrel        1.38     1.00       0.19

mod.coef.scaled %>% dplyr::group_by(sp) %>%
  reframe(range.slope=range(Estimate), range.r2=range(r2), range.rmse=range(rmse)) %>% as.data.frame() %>% mutate_if(is.numeric, round, digits=2)

#               sp range.slope range.r2 range.rmse
# 5          mouse        0.76     0.50       0.00
# 6          mouse        1.00     1.00       0.61

# 1       chipmunk        0.98     0.95       0.00
# 2       chipmunk        1.00     1.00       0.19

# 3 flyingsquirrel        0.95     0.90       0.00
# 4 flyingsquirrel        1.00     1.00       0.28

###############
#skip to here
###############
#final plots
allscr.scaled <- read.table(file="Data_linearregs/CR_estimates.txt", sep=",", header=T) #point estimates per grid per model per species
# allscr.scaled <- read.table(file="Data_linearregs/CR_estimates_correctgridnums.txt", sep=",", header=T) #point estimates per grid per model per species

mod.coef.all <- read.table(file="Data_linearregs/CR_slopes_individualidentity.txt", sep=",", header=T)
mod.coef.all.sc <- read.table(file="Data_linearregs/CR_slopes_spatialcount.txt", sep=",", header=T)
mod.coef.scaled <- read.table(file="Data_linearregs/CR_slopes_scaled.txt", sep=",", header=T)

##############

mod.coef.all$mod <- factor(mod.coef.all$mod, levels=order.mods.names)
mod.coef.all$sp <- factor(mod.coef.all$sp, levels=order.sp.names)

mod.coef.all.sc$mod <- factor(mod.coef.all.sc$mod, levels=order.mods.names)
mod.coef.all.sc$sp <- factor(mod.coef.all.sc$sp, levels=order.sp.names)

allscr.scaled$mod <- factor(allscr.scaled$mod, levels=order.mods.names, labels=order.mods.names)
allscr.scaled$model <- factor(allscr.scaled$model, levels=order.mods.names)
allscr.scaled$sp <- factor(allscr.scaled$sp, levels=order.sp.names, labels=order.sp.names)

# mod.coef.scaled.sc <- mod.coef.scaled[mod.coef.scaled$mod %in% c("SC info pool","SC no.info","SC info"),] #spatial count slopes
mod.coef.scaled <- mod.coef.scaled[!mod.coef.scaled$mod %in% c("SC info pool","SC no.info","SC info"),] #non spatial count slopes

mod.coef.scaled$mod <- factor(mod.coef.scaled$mod, levels=order.mods.names)
mod.coef.scaled$sp <- factor(mod.coef.scaled$sp, levels=order.sp.names)

# mod.coef.scaled.sc[mod.coef.scaled.sc$sp == "flyingsquirrel",]$sp <- "flying\nsquirrel"
# mod.coef.scaled.sc$sp <- factor(mod.coef.scaled.sc$sp, levels=order.sp.names)
# mod.coef.scaled.sc$mod <- factor(mod.coef.scaled.sc$mod, levels=order.mods.names)

all <-  ggplot(data=allscr.scaled[!allscr.scaled$model %in% c("SC no.info","SC info","SC info pool"),] , aes(x=dhat, y=dhat.random, col=model)) +
  geom_abline(slope=1, intercept=0, col="grey60", linewidth=0.25) +
  geom_errorbar(aes(ymin=X2.5..random, ymax=X97.5..random), width=0.01, lty="dashed", linewidth=0.1) + #vertical error bars
  geom_errorbar(aes(xmin=X2.5., xmax=X97.5.), width=0.01, lty="dashed", linewidth=0.1) + #horizontal error bars
  geom_point(size=3, pch=1, stroke=2) + #, alpha=0.5) +
  geom_smooth(method="lm", aes(fill=model, lty=model), alpha=0.1) + #add ribbons
  scale_color_manual(values=c.vals) + scale_fill_manual(values=c.vals) +
  # scale_x_log10() + scale_y_log10() +
  # scale_x_sqrt() + scale_y_sqrt() +
  xlab("model estimate (animals/ha)") + ylab("SCR random (animals/ha)") +
  theme_bw(base_size=20) + theme(panel.grid = element_blank(), legend.position="top", strip.background = element_blank(), strip.text.x = element_blank()) + facet_wrap(~sp, scales="free")
all

#slope
s <- ggplot(mod.coef.all, aes(x=mod, y=Estimate, col=mod)) + 
  geom_point(size=5) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1) +
  scale_color_manual(values=c.vals) +
  geom_hline(yintercept = 0, lty="dashed") +
  geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("slope") + scale_y_continuous(limits=c(0,NA)) +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + facet_wrap(~sp, scale="free_y")

#r2
r2 <- ggplot(mod.coef.all, aes(x=mod, y=r2, col=mod)) + geom_point(size=5) + 
  scale_color_manual(values=c.vals) +
  geom_hline(yintercept = 0, lty="dashed") + geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("R2") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + facet_wrap(~sp, scale="free_y")

#rmse
rmse <- ggplot(mod.coef.all, aes(x=mod, y=rmse, col=mod)) + geom_point(size=5) + 
  scale_color_manual(values=c.vals) +
  geom_hline(yintercept = 0, lty="dashed") + #geom_hline(yintercept=1, lty="dashed") + 
  xlab("model") + ylab("RMSE") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5)) + facet_wrap(~sp, scale="free_y")

p <- ggarrange(all, s, r2, rmse, ncol=1, align="v", heights=c(0.5,0.25,0.15,0.25))
p
# ggsave(p, filename="Figures/CR_fig.tiff", height=15, width=15, units="in", dpi=400, compression="lzw")

############
# all.sc <-  ggplot(data=allscr.scaled[allscr.scaled$model %in% c("SC no.info","SC info","SC info pool"),] , aes(x=dhat, y=dhat.random, col=model)) +
#   geom_abline(slope=1, intercept=0, col="grey60", linewidth=0.25) +
#   geom_errorbar(aes(ymin=X2.5..random, ymax=X97.5..random), width=0.01, lty="dashed", linewidth=0.1) + #vertical error bars
#   geom_errorbar(aes(xmin=X2.5., xmax=X97.5.), width=0.01, lty="dashed", linewidth=0.1) + #horizontal error bars
#   geom_point(size=3, pch=1, stroke=2) + #, alpha=0.5) +
#   geom_smooth(method="lm", aes(fill=model, lty=model), alpha=0.1) + #add ribbons
#   scale_color_manual(values=sc.vals) + scale_fill_manual(values=sc.vals) +
#   # scale_x_log10() + scale_y_log10() +
#   # scale_x_sqrt() + scale_y_sqrt() +
#   xlab("model estimate (animals/ha)") + ylab("SCR random (animals/ha)") +
#   theme_bw(base_size=20) + theme(panel.grid = element_blank(), legend.position="top", strip.background = element_blank(), strip.text.x = element_blank()) + facet_wrap(~sp, scales="free")
# all.sc
# 
# #slope
# s.sc <- ggplot(mod.coef.all.sc, aes(x=mod, y=Estimate, col=mod)) + 
#   geom_point(size=5) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1) +
#   scale_color_manual(values=sc.vals) +
#   geom_hline(yintercept = 0, lty="dashed") +
#   geom_hline(yintercept=1, lty="dashed") + 
#   xlab("") + ylab("slope") + scale_y_continuous(limits=c(0,NA)) +
#   theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + facet_wrap(~sp, scale="free_y")
# 
# #r2
# r2.sc <- ggplot(mod.coef.all.sc, aes(x=mod, y=r2, col=mod)) + geom_point(size=5) + 
#   scale_color_manual(values=sc.vals) +
#   geom_hline(yintercept = 0, lty="dashed") + geom_hline(yintercept=1, lty="dashed") + 
#   xlab("") + ylab("R2") +
#   theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + facet_wrap(~sp, scale="free_y")
# 
# #rmse
# rmse.sc <- ggplot(mod.coef.all.sc, aes(x=mod, y=rmse, col=mod)) + geom_point(size=5) + 
#   scale_color_manual(values=sc.vals) +
#   geom_hline(yintercept = 0, lty="dashed") + #geom_hline(yintercept=1, lty="dashed") + 
#   xlab("model") + ylab("RMSE") +
#   theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5)) + facet_wrap(~sp, scale="free_y")
# 
# p.sc <- ggarrange(all.sc, s.sc, r2.sc, rmse.sc, ncol=1, align="v", heights=c(0.5,0.25,0.15,0.25))
# p.sc
# # ggsave(p.sc, filename="Figures/CR_sc_fig.tiff", height=15, width=15, units="in", dpi=400, compression="lzw")

##################################
#for scaled values
##################################

# all.s <-  ggplot(data=allscr.scaled, aes(x=dhat.scaled, y=dhat.random.scaled, col=model)) +
all.s <-  ggplot(data=allscr.scaled[!allscr.scaled$mod %in% c("SC info pool","SC no.info","SC info"),], aes(x=dhat.scaled, y=dhat.random, col=model)) +
  # geom_abline(slope=1, intercept=0, col="grey60", linewidth=0.25) +
  # geom_line(data=allscr.scaled[allscr.scaled$model == "SCR random",], aes(x=dhat.scaled, y=dhat.random), col="grey60", linewidth=0.25) +
  geom_smooth(method="lm", data=allscr.scaled[allscr.scaled$model == "SCR random",], aes(x=dhat.scaled, y=dhat.random), linewidth=0.25, col="grey80", fullrange=T) +
  geom_errorbar(aes(ymin=X2.5..random, ymax=X97.5..random), width=0.01, lty="dashed", linewidth=0.1) + #vertical error bars
  geom_errorbar(aes(xmin=X2.5.scaled, xmax=X97.5.scaled), width=0.01, lty="dashed", linewidth=0.1) + #horizontal error bars
  geom_point(size=3, pch=1, stroke=2) + #, alpha=0.5) +
  geom_smooth(method="lm", aes(fill=model, lty=model), alpha=0.1) + #add ribbons
  scale_color_manual(values=c.vals) + scale_fill_manual(values=c.vals) +
  # scale_x_log10() + scale_y_log10() +
  # scale_x_sqrt() + scale_y_sqrt() +
  xlab("scaled model estimate (animals/ha)") + ylab("SCR random (animals/ha)") +
  theme_bw(base_size=20) + theme(panel.grid = element_blank(), legend.position="top", strip.background = element_blank(), strip.text.x = element_blank()) + facet_wrap(~sp, scales="free")
all.s

#slope
s.s <- ggplot(mod.coef.scaled, aes(x=mod, y=Estimate, col=mod)) + 
  geom_point(size=5) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1) +
  scale_color_manual(values=c.vals) +
  geom_hline(yintercept = 0, lty="dashed") +
  geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("slope") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + facet_wrap(~sp, scale="free_y")

#r2
r2.s <- ggplot(mod.coef.scaled, aes(x=mod, y=r2, col=mod)) + geom_point(size=5) + 
  scale_color_manual(values=c.vals) +
  geom_hline(yintercept = 0, lty="dashed") + geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("R2") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + facet_wrap(~sp, scale="free_y")

#rmse
rmse.s <- ggplot(mod.coef.scaled, aes(x=mod, y=rmse, col=mod)) + geom_point(size=5) + 
  scale_color_manual(values=c.vals) +
  geom_hline(yintercept = 0, lty="dashed") + #geom_hline(yintercept=1, lty="dashed") + 
  xlab("model") + ylab("RMSE") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5)) + facet_wrap(~sp, scale="free_y")

p.s <- ggarrange(all.s, s.s, r2.s, rmse.s, ncol=1, align="v", heights=c(0.5,0.25,0.15,0.25))
p.s
ggsave(p.s, filename="Figures/CR_scaled_fig.tiff", height=15, width=15, units="in", dpi=400, compression="lzw")

###################
#mouse, chip, flyer
#spatial count models with Capture recapture data
allscr.scaled$model <- factor(allscr.scaled$model, levels=order.mods.names)

cr.sc.s <-  ggplot(data=allscr.scaled[allscr.scaled$mod %in% c("SC info pool","SC no.info","SC info"),], aes(x=dhat.scaled, y=dhat.random, col=model)) +
  # geom_abline(slope=1, intercept=0, col="grey60", linewidth=0.25) +
  # geom_line(data=allscr.scaled[allscr.scaled$model == "SCR random",], aes(x=dhat.scaled, y=dhat.random), col="grey60", linewidth=0.25) +
  geom_smooth(method="lm", data=allscr.scaled[allscr.scaled$model == "SCR random",], aes(x=dhat.scaled, y=dhat.random), linewidth=0.25, col="grey80", fullrange=T) +
  geom_errorbar(aes(ymin=X2.5..random, ymax=X97.5..random), width=0.01, lty="dashed", linewidth=0.1) + #vertical error bars
  # geom_errorbar(aes(xmin=X2.5.scaled, xmax=X97.5.scaled), width=0.01, lty="dashed", linewidth=0.1) + #horizontal error bars
  geom_point(size=3, pch=1, stroke=2) + #, alpha=0.5) +
  geom_smooth(method="lm", aes(fill=model, lty=model), alpha=0.1) + #add ribbons
  scale_color_manual(values=sc.vals) + scale_fill_manual(values=sc.vals) +
  # scale_x_log10() + scale_y_log10() +
  # scale_x_sqrt() + scale_y_sqrt() +
  xlab("scaled model estimate (animals/ha)") + ylab("SCR random (animals/ha)") + #coord_cartesian(xlim=c(NA, 5)) +
  theme_bw(base_size=20) + theme(panel.grid = element_blank(), legend.position="top", strip.background = element_blank(), strip.text.x = element_blank()) + 
  facet_wrap(~sp, scales="free")
hbars <- allscr.scaled[allscr.scaled$mod %in% c("SC info pool","SC no.info","SC info"),] #because error bars make chipmunk and flyer graphs are weird
hbars$mod <- factor(hbars$mod, levels=order.mods.names)
hbars[hbars$sp == "chipmunk" & hbars$X97.5.scaled > 5,]$X97.5.scaled <- 5
hbars[hbars$sp == "flying\nsquirrel" & hbars$X97.5.scaled > 5,]$X97.5.scaled <- 5
cr.sc.s <- cr.sc.s + geom_errorbar(data=hbars, aes(xmin=X2.5.scaled, xmax=X97.5.scaled), width=0.01, lty="dashed", linewidth=0.1) #horizontal error bars
cr.sc.s

#slope
s.sc.s <- ggplot(mod.coef.scaled.sc, aes(x=mod, y=Estimate, col=mod)) + 
  geom_point(size=5) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1) +
  scale_color_manual(values=sc.vals) +
  geom_hline(yintercept = 0, lty="dashed") +
  geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("slope") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + 
  facet_wrap(~sp, scale="free_y")

#r2
r2.sc.s <- ggplot(mod.coef.scaled.sc, aes(x=mod, y=r2, col=mod)) + geom_point(size=5) + 
  scale_color_manual(values=sc.vals) +
  geom_hline(yintercept = 0, lty="dashed") + geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("R2") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + facet_wrap(~sp, scale="free_y")

#rmse
rmse.sc.s <- ggplot(mod.coef.scaled.sc, aes(x=mod, y=rmse, col=mod)) + geom_point(size=5) + 
  scale_color_manual(values=sc.vals) +
  geom_hline(yintercept = 0, lty="dashed") + #geom_hline(yintercept=1, lty="dashed") + 
  xlab("model") + ylab("RMSE") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5)) + facet_wrap(~sp, scale="free_y")

p.sc.s <- ggarrange(cr.sc.s, s.sc.s, r2.sc.s, rmse.sc.s, ncol=1, align="v", heights=c(0.5,0.25,0.15,0.25))
p.sc.s
ggsave(p.sc.s, filename="Figures/CR_sc_scaled_fig.tiff", height=15, width=15, units="in", dpi=400, compression="lzw")

################################
#summaries of each model

table(allscr.scaled[,c("sp","model")])
#                    model
# sp                 MNKA Huggins SCR separate SCR pool SCR random SC no.info SC info SC info pool
# mouse               8       7            8        8          8          0       8            8
# chipmunk            9       9            9        9          9          1       9            9
# flying\nsquirrel    9       9            9        9          9          2       9            9

#missing grid 6 and grid 9 for Huggins PEMA

library(dplyr)
mod.coef.scaled %>% 
  mutate_if(is.numeric, round, digits=3)

#    Estimate Std..Error      t.value Pr...t..             sp          mod    r2  rmse
# 1     1.000      0.000 2.061796e+16    0.000          mouse   SCR random 1.000 0.000
# 2     0.909      0.170 5.349000e+00    0.002          mouse     SCR pool 0.798 0.389
# 3     0.982      0.077 1.275500e+01    0.000          mouse SCR separate 0.959 0.176
# 4     0.834      0.247 3.376000e+00    0.020          mouse      Huggins 0.634 0.511
# 5     0.990      0.058 1.704100e+01    0.000          mouse         MNKA 0.976 0.133

# 8     1.000      0.000 1.577801e+16    0.000 flyingsquirrel   SCR random 1.000 0.000
# 9     0.975      0.084 1.160300e+01    0.000 flyingsquirrel     SCR pool 0.944 0.210
# 10    0.979      0.078 1.262700e+01    0.000 flyingsquirrel SCR separate 0.952 0.193
# 11    0.970      0.092 1.056000e+01    0.000 flyingsquirrel      Huggins 0.932 0.229
# 12    0.919      0.149 6.154000e+00    0.000 flyingsquirrel         MNKA 0.822 0.372

# 16    1.000      0.000 1.230297e+16    0.000       chipmunk   SCR random 1.000 0.000
# 17    0.985      0.065 1.515200e+01    0.000       chipmunk     SCR pool 0.966 0.162
# 18    0.998      0.021 4.754400e+01    0.000       chipmunk SCR separate 0.996 0.052
# 19    0.980      0.075 1.299800e+01    0.000       chipmunk      Huggins 0.955 0.188
# 20    0.874      0.183 4.768000e+00    0.002       chipmunk         MNKA 0.731 0.457

mod.coef.scaled.sc %>%
  mutate_if(is.numeric, round, digits=3)
#    Estimate Std..Error t.value Pr...t..               sp          mod     r2  rmse
# 6     0.621      0.320   1.939    0.110            mouse SC info pool  0.315 0.699
# 7     0.811      0.354   2.288    0.071            mouse      SC info  0.414 0.647

# 13    0.258      0.365   0.705    0.503 flying\nsquirrel SC info pool -0.067 0.911
# 14    1.000         NA      NA       NA flying\nsquirrel   SC no.info     NA 0.000
# 15    0.665      0.282   2.353    0.051 flying\nsquirrel      SC info  0.362 0.705

# 21    0.718      0.263   2.728    0.029         chipmunk SC info pool  0.446 0.656
# 22    0.794      0.230   3.457    0.011         chipmunk      SC info  0.578 0.573