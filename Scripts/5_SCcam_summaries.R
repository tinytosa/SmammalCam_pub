#########################
#########################

#SCR random on the y axis
#SC cam estimates on x axis

#########################
#load packages
require(dplyr)
require(plyr) # for ddply
require(tidyr)

require(ggplot2)
require(ggpubr)
require(viridis)

#########################
#load data

#conversion table
sp.convert <- data.frame(sp=c("chipmunk","flying\nsquirrel","flying\nsquirrel","mouse"), acro=c("TATO","GLSA","GLOR","PEMA"))
# sp.convert <- data.frame(sp=c("chipmunk","flying\nsquirrel","mouse"), sp2=c("TOWNSENDS CHIPMUNK","FLYING SQUIRREL","DEER MOUSE"))
sp.convert$sp <- factor(sp.convert$sp, levels=c("mouse","chipmunk","flying\nsquirrel"))

order.mods <- c("CAM-SC_noinfo_separate", "CAM-SC_info_separate", "CAM-SC_noinfo_multigrid", "CAM-SC_info_multigrid")
order.mods.names <- c("SC no.info sep","SC info sep", "SC no.info pool", "SC info pool")

#random SCR estimates
# random <- read.table(file="Data_linearregs/CR_estimates.txt", sep=",", header=T) #wrong grid numbers for PEMA
random <- read.table(file="Data_linearregs/CR_estimates_correctgridnums.txt", sep=",", header=T)
random <- random[random$model == "SCR random",] #has absolute values and scaled values
random <- random[,c("sp","grid",
                    "dhat","mean.random","sd.random","X2.5..random","X97.5..random", #remove extra columns
                    "dhat.scaled","X2.5.r.scaled","X97.5.r.scaled")]
names(random) <- c("sp","grid", "dhat.random","mean.random","sd.random","X2.5..random","X97.5..random","dhat.random.scaled","X2.5.r.scaled","X97.5.r.scaled")
table(random$sp)

random[random$sp == "flyingsquirrel",]$sp <- "flying\nsquirrel"

#
# sc.cam <- read.table(file="Output_cam/output_CAM-SC_estimates_2023-06-13.txt", sep=",", header=T)
sc.cam <- read.table(file="Output_cam/output_CAM-SC_separate_estimates.txt", sep=",", header=T) #separate
sc.cam2 <- read.table(file="Output_cam/output_CAM-SC_multigrid_estimates.txt", sep=",", header=T) #multigrid

#match sc.cam2 columns to sc.cam
sc.cam2$model <- sc.cam2$model.y
sc.cam2$mode <- sc.cam2$est


sc.cam <- rbind(sc.cam, sc.cam2[,names(sc.cam2) %in% names(sc.cam)])
sc.cam <- sc.cam[grep(sc.cam$param, pattern="D"),]

table(sc.cam[,c("t","sp","model")])

sc.cam <- sc.cam[!duplicated(sc.cam),] #check for and remove duplicates
sc.cam <- sc.cam[!is.na(sc.cam$Rhat),] #remove rows with NA for rhat
sc.cam <- sc.cam[sc.cam$Rhat < 1.2,] #remove rows with high rhat values
sc.cam$t <- factor(sc.cam$t, levels=c(0,15,60,1440))

table(sc.cam[,c("t","sp","model")])
# , , model = CAM-SC_info_multigrid
# 
# sp
# t      GLSA PEMA TATO
# 0       7    6    8
# 15      0    6    0
# 60      0    8    0
# 1440    0    6    7
# 
# , , model = CAM-SC_info_separate
# 
# sp
# t      GLSA PEMA TATO
# 0       2    0    1
# 15      4    1    8
# 60      4    3    6
# 1440    7    2    5
# 
# , , model = CAM-SC_noinfo_multigrid
# 
# sp
# t      GLSA PEMA TATO
# 0       0    8    3
# 15      3    4    8
# 60      2    7    8
# 1440    0    0    0
# 
# , , model = CAM-SC_noinfo_separate
# 
# sp
# t      GLSA PEMA TATO
# 0       3    0    2
# 15      3    0    5
# 60      3    0    1
# 1440    2    0    0

#######################################

sc.cam <- merge(sp.convert, sc.cam, by.y="sp", by.x="acro", all.y=T)
sc.cam <- merge(sc.cam, random, by.x=c("g","sp"), by.y=c("grid","sp"), all.x=T)

#need to remove grid 6 PEMA from analysis since no capture-recapture
sc.cam <- sc.cam[!(sc.cam$g == 6 & sc.cam$sp == "mouse"),]

sc.cam$t <- factor(sc.cam$t, levels=c(0,15,60,1440))
# pivot_wider(plyr::ddply(sc.cam, .(sp, t, model), summarize, count=length(g)), values_from="count", names_from="sp")

table(sc.cam[,c("t","sp","model")])
# , , model = CAM-SC_info_multigrid
# 
#       sp
# t      mouse chipmunk flying\nsquirrel
# 0        5        8                7
# 15       5        0                0
# 60       7        0                0
# 1440     6        7                0
# 
# , , model = CAM-SC_info_separate
# 
# sp
# t      mouse chipmunk flying\nsquirrel
# 0        0        1                2
# 15       1        8                4
# 60       3        6                4
# 1440     2        5                7
# 
# , , model = CAM-SC_noinfo_multigrid
# 
# sp
# t      mouse chipmunk flying\nsquirrel
# 0        7        3                0
# 15       3        8                3
# 60       6        8                2
# 1440     0        0                0
# 
# , , model = CAM-SC_noinfo_separate
# 
# sp
# t      mouse chipmunk flying\nsquirrel
# 0        0        2                3
# 15       0        5                3
# 60       0        1                3
# 1440     0        0                2

#scale sc.cam mean values
#only scale values that converged
#create new columns
sc.cam$mode.scaled <- NA
sc.cam$X2.5.scaled <- NA
sc.cam$X97.5.scaled <- NA

mod.coef.all <- data.frame()
mod.coef.scaled <- data.frame()

for(t in unique(sc.cam$t))
{
  print(t)
  for(sp in unique(sc.cam$sp))
  {
    print(sp)
    for(m in unique(sc.cam$model))
    {
      print(m)
      if(nrow(sc.cam[sc.cam$t == t & sc.cam$sp == sp & sc.cam$model == m,]) > 0)
      {
        ########
        #not scaled
        #linear regression with non-scaled values
        l <- lm(dhat.random ~ mode, data=sc.cam[sc.cam$t == t & sc.cam$sp == sp & sc.cam$model == m,])
        
        mod.coef <- data.frame(coef(summary(l)))[2,]
        mod.coef$sp <- sp
        mod.coef$t <- t
        mod.coef$mod <- m
        mod.coef$r2 <- summary(l)$adj.r.squared
        mod.coef$rmse <- sqrt(mean(l$residuals^2))
        mod.coef$n <- nrow(sc.cam[sc.cam$t == t & sc.cam$sp == sp & sc.cam$model == m,])
        mod.coef.all <- rbind(mod.coef.all, mod.coef)
        
        ########
        #scaled
        s <- scale(sc.cam[sc.cam$t == t & sc.cam$sp == sp  & sc.cam$model == m,]$mode)
        sc.cam[sc.cam$t == t & sc.cam$sp == sp & sc.cam$model == m,]$mode.scaled <- s
        
        sc.cam[sc.cam$t == t & sc.cam$sp == sp & sc.cam$model == m,]$X2.5.scaled <- (sc.cam[sc.cam$t == t & sc.cam$sp == sp & sc.cam$model == m,]$X2.5. - attr(s, "scaled:center"))/attr(s, "scaled:scale")
        sc.cam[sc.cam$t == t & sc.cam$sp == sp & sc.cam$model == m,]$X97.5.scaled <- (sc.cam[sc.cam$t == t & sc.cam$sp == sp & sc.cam$model == m,]$X97.5. - attr(s, "scaled:center"))/attr(s, "scaled:scale")
        
        if(nrow(sc.cam[sc.cam$t == t & sc.cam$sp == sp & sc.cam$model == m,]) == 1)
        {
          sc.cam[sc.cam$t == t & sc.cam$sp == sp & sc.cam$model == m,]$mode.scaled <- 0
        }
        
        l <- lm(dhat.random.scaled ~ mode.scaled, data=sc.cam[sc.cam$t == t & sc.cam$sp == sp & sc.cam$model == m,]) #linear regression
        
        mod.coef.s <- data.frame(coef(summary(l)))[2,]
        mod.coef.s$sp <- sp
        mod.coef.s$t <- t
        mod.coef.s$mod <- m
        mod.coef.s$r2 <- summary(l)$adj.r.squared
        mod.coef.s$rmse <- sqrt(mean(l$residuals^2))
        mod.coef.s$n <- nrow(sc.cam[sc.cam$t == t & sc.cam$sp == sp & sc.cam$model == m,])
        mod.coef.scaled <- rbind(mod.coef.scaled, mod.coef.s)
      }
    }
  }
}

# write.table(sc.cam, file="Data_linearregs/sc.cam_estimates.txt", sep=",", row.names=F)

#################
#1b. plot unscaled values together
random %>% group_by(sp) %>% summarize(min(dhat.random), max(dhat.random), min(X2.5..random), max(X97.5..random)) #summarize min and max values of SCR random
empty <- data.frame(sp=c("mouse","mouse","chipmunk","chipmunk","flying\nsquirrel","flying\nsquirrel"), mode=rep(0,6), dhat.random=c(9.16, 96.0, 1.07, 10.4, 0.0662, 4.28), model=rep("",6), t=rep(0,6))
empty$sp <- factor(empty$sp, levels=c("mouse","chipmunk","flyingsquirrel","flying\nsquirrel"))
empty$t <- factor(empty$t, levels=c(0,15,60,1440))

mod.coef.all$mod <- factor(mod.coef.all$mod, levels=order.mods, labels=order.mods.names)
mod.coef.all$sp <- factor(mod.coef.all$sp, levels=c("mouse","chipmunk","flyingsquirrel","flying\nsquirrel"))
mod.coef.all$t <- factor(mod.coef.all$t, levels=c(0,15,60,1440))

sc.cam$model <- factor(sc.cam$model, levels=order.mods, labels=order.mods.names)
sc.cam$t <- factor(sc.cam$t, levels=c(0,15,60,1440))

# all.o <- ggplot(data=empty, aes(x=mode, y=dhat.random, col=t, fill=t, shape=model)) +
#   geom_point(pch=NA) +
all.o <- ggplot(data=sc.cam, aes(x=mode, y=dhat.random, col=t, fill=t, shape=model)) +
  scale_shape_manual(name="model", values=c(2, 17, 1, 16)) +
  # geom_smooth(method="lm", aes(lty=model), alpha=0.1) + #add ribbons
  # geom_errorbar(aes(ymin=X2.5..random, ymax=X97.5..random), width=0.01, lty="dashed", linewidth=0.1) + #vertical error bars
  # geom_errorbar(aes(xmin=X2.5., xmax=X97.5.), width=0.01, lty="dashed", linewidth=0.1) +
  geom_point(size=3, stroke=2) + #, alpha=0.5) +
  # scale_color_manual(values=sc.vals) + scale_fill_manual(values=sc.vals) +
  # scale_linetype(name="model") +
  # scale_shape_manual(name="model", values=c(2, 17, 1, 16)) +
  facet_wrap(~sp, scales="free") +
  theme_bw(base_size=20) + theme(panel.grid = element_blank(), legend.position="top", strip.background = element_blank(), strip.text.x = element_blank()) +
  coord_cartesian(ylim=c(0, NA))
all.o

all.0 <- ggplot(data=sc.cam, aes(x=mode, y=dhat.random, col=t, fill=t, shape=model)) + 
  geom_abline(aes(slope=1, intercept=0), col="grey50", lwd=0.5) + 
  geom_smooth(method="lm", aes(lty=model), alpha=0.1, se=F) + #add ribbons
  geom_errorbar(aes(ymin=X2.5..random, ymax=X97.5..random), width=0.01, lty="dashed", linewidth=0.1) + #vertical error bars
  geom_errorbar(aes(xmin=X2.5., xmax=X97.5.), width=0.01, lty="dashed", linewidth=0.1) +
  geom_point(size=3, stroke=2) + 
  scale_shape_manual(name="model", values=c(2, 17, 1, 16)) +
  scale_linetype(name="model") +
  # facet_wrap(~t + sp, scales="free", ncol=3) +
  facet_wrap(~sp, scales="free", ncol=3) +
  xlab("CAM SC density estimate (animals/ha)") + ylab("SCR random (animals/ha)") +
  theme_bw(base_size=20) + theme(panel.grid = element_blank(), legend.position="top", strip.background = element_blank(), strip.text.x = element_blank()) +
  coord_cartesian(ylim=c(0, NA))
all.0
# ggsave(all.0, file="Figures/CAM-SC_1.tiff", height=20, width=12, units="in", dpi=400, compression="lzw")
# ggsave(all.0, file="Figures/CAM-SC_2.tiff", height=6, width=12, units="in", dpi=400, compression="lzw")

#slope
s <- ggplot(mod.coef.all, aes(x=mod, y=Estimate, col=t, shape=mod)) + 
  geom_point(size=5, position=position_dodge(width=1), stroke=2) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=1)) +
  # scale_color_manual(values=sc.vals) +
  scale_shape_manual(name="model", values=c(2, 17, 1, 16)) +
  geom_hline(yintercept=0, lty="dashed") +
  geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("slope") + #scale_y_continuous(limits=c(0,NA)) +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + 
  facet_wrap(~sp, scale="free_y")

#r2
r2 <- ggplot(mod.coef.all, aes(x=mod, y=r2, col=t, shape=mod)) + geom_point(size=5, position=position_dodge(width=1), stroke=2) + 
  # scale_color_manual(values=sc.vals) +
  scale_shape_manual(name="model", values=c(2, 17, 1, 16)) +
  geom_hline(yintercept = 0, lty="dashed") + geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("R2") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + 
  facet_wrap(~sp, scale="free_y")

#rmse
rmse <- ggplot(mod.coef.all, aes(x=mod, y=rmse, col=t, shape=mod)) + geom_point(size=5, position=position_dodge(width=1), stroke=2) + 
  # scale_color_manual(values=sc.vals) +
  scale_shape_manual(name="model", values=c(2, 17, 1, 16)) +
  geom_hline(yintercept = 0, lty="dashed") + #geom_hline(yintercept=1, lty="dashed") + 
  xlab("model") + ylab("RMSE") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5)) + facet_wrap(~sp, scale="free_y")

p <- ggarrange(all.0, s, r2, rmse, ncol=1, align="v", heights=c(0.5,0.25,0.15,0.25))
p
# ggsave(p, filename="Figures/CAM-SC_fig.tiff", height=15, width=15, units="in", dpi=400, compression="lzw")

################
#1c. scale values and plot and calculate slopes, RMSE, and R2
mod.coef.scaled$sp <- factor(mod.coef.scaled$sp, levels=c("mouse","chipmunk","flyingsquirrel","flying\nsquirrel"))
mod.coef.scaled$mod <- factor(mod.coef.scaled$mod, levels=order.mods, labels=order.mods.names)
mod.coef.scaled$t <- factor(mod.coef.scaled$t, levels=c(0,15,60,1440))

#create empty data frame with all possible combinations of species, t, model
e <- data.frame(sp=rep(c("mouse","chipmunk","flying\nsquirrel"), each=16),
                t=rep(rep(c(0,15,60,1440), each=4), times=3),
                mod=rep(c("SC no.info sep","SC info sep","SC no.info pool","SC info pool"), times=12))
e$sp <- factor(e$sp, levels=c("mouse","chipmunk","flyingsquirrel","flying\nsquirrel"))
e$t <- factor(e$t, levels=c(0,15,60,1440))
e$mod <- factor(e$mod, levels=order.mods.names)

m2 <- mod.coef.scaled
m <- merge(mod.coef.scaled, e, by=c("sp","t","mod"), all.y=T)
mod.coef.scaled <- m

all.0 <- ggplot(data=sc.cam, aes(x=mode.scaled, y=dhat.random, col=t, shape=model)) +
  # geom_smooth(data=sc.cam, method="lm", aes(fill=t, lty=model), alpha=0.1) + #add ribbons
  geom_smooth(data=sc.cam, method="lm", aes(lty=model), alpha=0.1, se=F) + #add ribbons
  geom_errorbar(aes(ymin=X2.5..random, ymax=X97.5..random), width=0.01, lty="dashed", linewidth=0.1) + #vertical error bars
  # geom_errorbar(aes(xmin=X2.5.scaled, xmax=X97.5.scaled), width=0.01, lty="dashed", linewidth=0.1) +
  geom_point(size=3, stroke=1) + #, alpha=0.5) +
  # scale_color_manual(values=sc.vals) + scale_fill_manual(values=sc.vals) +
  scale_shape_manual(name="model", values=c(2, 17, 1, 16)) +
  scale_linetype(name="model") +
  # scale_shape_manual(name="model", values=c(2, rep(1, 5))) +
  xlab("scaled model estimate (animals/ha)") + ylab("SCR random (animals/ha)") +
  theme_bw(base_size=20) + theme(panel.grid = element_blank(), legend.position="top", strip.background = element_blank(), strip.text.x = element_blank()) + 
  facet_wrap(~sp, scales="free") +
  coord_cartesian(ylim=c(0, NA))
# coord_cartesian(ylim=c(0, 150))
all.0

#slope
s <- ggplot(mod.coef.scaled, aes(x=mod, y=Estimate, col=t, shape=mod)) +
# s <- ggplot(mod.coef.scaled, aes(x=t, y=Estimate, col=t, shape=mod)) + 
  geom_point(size=5, position=position_dodge(width=0.75), stroke=2) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=0.75)) +
  # scale_color_manual(values=sc.vals) +
  scale_shape_manual(name="model", values=c(2, 17, 1, 16)) +
  geom_hline(yintercept=0, lty="dashed") +
  geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("slope") + #scale_y_continuous(limits=c(0,NA)) +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + 
  facet_wrap(~sp, scale="free_y")

#r2
r2 <- ggplot(mod.coef.scaled[mod.coef.scaled$n > 1 & !is.na(mod.coef.scaled$n),], aes(x=mod, y=r2, col=t, shape=mod)) + geom_point(size=5, position=position_dodge(width=0.75), stroke=2) +
# r2 <- ggplot(mod.coef.scaled[mod.coef.scaled$n > 1 & !is.na(mod.coef.scaled$n),], aes(x=t, y=r2, col=t, shape=mod)) + geom_point(size=5, position=position_dodge(width=0.75), stroke=2) + 
  # scale_color_manual(values=sc.vals) +
  scale_shape_manual(name="model", values=c(2, 17, 1, 16)) +
  geom_hline(yintercept = 0, lty="dashed") + geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("R2") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + 
  facet_wrap(~sp, scale="free_y")

#rmse
rmse <- ggplot(mod.coef.scaled[mod.coef.scaled$n > 1 & !is.na(mod.coef.scaled$n),], aes(x=mod, y=rmse, col=t, shape=mod)) + geom_point(size=5, position=position_dodge(width=0.75), stroke=2) +
# rmse <- ggplot(mod.coef.scaled[mod.coef.scaled$n > 1 & !is.na(mod.coef.scaled$n),], aes(x=t, y=rmse, col=t, shape=mod)) + geom_point(size=5, position=position_dodge(width=0.75), stroke=2) + 
  # scale_color_manual(values=sc.vals) +
  scale_shape_manual(name="model", values=c(2, 17, 1, 16)) +
  geom_hline(yintercept = 0, lty="dashed") + #geom_hline(yintercept=1, lty="dashed") + 
  xlab("model") + ylab("RMSE") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5)) + facet_wrap(~sp, scale="free_y")

p <- ggarrange(all.0, s, r2, rmse, ncol=1, align="v", heights=c(0.49,0.25,0.2,0.325))
p
# ggsave(p, filename="Figures/CAM-SC_scaled_fig.tiff", height=15, width=15, units="in", dpi=400, compression="lzw")

# write.csv(mod.coef.scaled, file="Data_linearregs/CAM-SC_slopes.csv")


##########
#create insets for slope for chipmunk
ggplot(mod.coef.scaled[mod.coef.scaled$sp == "chipmunk" & mod.coef.scaled$mod %in% c("SC no.info sep"),], aes(x=mod, y=Estimate, col=t, shape=mod)) +
  # s <- ggplot(mod.coef.scaled, aes(x=t, y=Estimate, col=t, shape=mod)) + 
  geom_point(size=5, position=position_dodge(width=0.75), stroke=2) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=0.75)) +
  # scale_color_manual(values=sc.vals) +
  scale_shape_manual(name="model", values=c(2, 17, 1, 16)) +
  geom_hline(yintercept=0, lty="dashed") +
  geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("") + #scale_y_continuous(limits=c(0,NA)) +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "null")) + 
  facet_wrap(~sp, scale="free_y")
ggsave(filename = "Figures/CAM-SC_scaled_fig_s1.tiff", height=2.5, width=3*3/5, unit="in", compression="lzw", dpi=400)

ggplot(mod.coef.scaled[mod.coef.scaled$sp == "chipmunk" & mod.coef.scaled$mod %in% c("SC no.info pool","SC info pool"),], aes(x=mod, y=Estimate, col=t, shape=mod)) +
  # s <- ggplot(mod.coef.scaled, aes(x=t, y=Estimate, col=t, shape=mod)) + 
  geom_point(size=5, position=position_dodge(width=0.75), stroke=2) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=0.75)) +
  # scale_color_manual(values=sc.vals) +
  scale_shape_manual(name="model", values=c(1, 16)) +
  geom_hline(yintercept=0, lty="dashed") +
  geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("") + scale_y_continuous(limits=c(0,NA), breaks=seq(0,1, 0.2)) +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "null")) + 
  facet_wrap(~sp, scale="free_y")
ggsave(filename = "Figures/CAM-SC_scaled_fig_s2.tiff", height=2.5, width=3.25, unit="in", compression="lzw", dpi=400)

ggplot(mod.coef.scaled[mod.coef.scaled$sp == "flying\nsquirrel" & mod.coef.scaled$mod %in% c("SC no.info sep"),], aes(x=mod, y=Estimate, col=t, shape=mod)) +
  # s <- ggplot(mod.coef.scaled, aes(x=t, y=Estimate, col=t, shape=mod)) + 
  geom_point(size=5, position=position_dodge(width=0.75), stroke=2) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=0.75)) +
  # scale_color_manual(values=sc.vals) +
  scale_shape_manual(name="model", values=c(2, 17, 1, 16)) +
  geom_hline(yintercept=0, lty="dashed") +
  geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("") + #scale_y_continuous(limits=c(0,NA)) +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "null")) + 
  facet_wrap(~sp, scale="free_y")
ggsave(filename = "Figures/CAM-SC_scaled_fig_s3.tiff", height=2.5, width=3*3/5, unit="in", compression="lzw", dpi=400)

ggplot(mod.coef.scaled[mod.coef.scaled$sp == "flying\nsquirrel" & mod.coef.scaled$mod %in% c("SC no.info pool","SC info pool"),], aes(x=mod, y=Estimate, col=t, shape=mod)) +
  # s <- ggplot(mod.coef.scaled, aes(x=t, y=Estimate, col=t, shape=mod)) + 
  geom_point(size=5, position=position_dodge(width=0.75), stroke=2) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=0.75)) +
  # scale_color_manual(values=sc.vals) +
  scale_shape_manual(name="model", values=c(1, 16)) +
  geom_hline(yintercept=0, lty="dashed") +
  geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("") + 
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "null")) + 
  facet_wrap(~sp, scale="free_y")
ggsave(filename = "Figures/CAM-SC_scaled_fig_s4.tiff", height=2.5, width=3.25, unit="in", compression="lzw", dpi=400)


# range can do the same thing as min and max without having to type it all out.
# mod.coef.scaled %>% group_by(sp) %>%
#   summarize(min(Estimate, na.rm=T), max(Estimate, na.rm=T), max(r2,na.rm=T), min(rmse,na.rm=T))

#ranges of slope, R2, and RMSE
detach(package:plyr)
mod.coef.scaled[!is.na(mod.coef.scaled$Estimate),] %>% dplyr::group_by(sp) %>%
  reframe(range(Estimate, na.rm=T), range(r2,na.rm=T), range(rmse,na.rm=T)) %>% as.data.frame() %>% mutate_if(is.numeric, round, digits=2)

#                 sp range(Estimate, na.rm = T) range(r2, na.rm = T) range(rmse, na.rm = T)
# 1            mouse                      -0.97                -0.19                   0.00
# 2            mouse                       1.30                 0.95                   0.96
# 3         chipmunk                      -0.84                -0.33                   0.00
# 4         chipmunk                       0.87                 1.00                   0.89
# 5 flying\nsquirrel                      -1.01                -0.48                   0.00
# 6 flying\nsquirrel                       0.88                 0.66                   0.79

mod.coef.scaled[mod.coef.scaled$rmse > 0 & !is.na(mod.coef.scaled$Estimate),] %>% dplyr::group_by(sp) %>%
  reframe(range(r2,na.rm=T), range(rmse,na.rm=T)) %>% as.data.frame() %>% mutate_if(is.numeric, round, digits=2)
#                 sp range(r2, na.rm = T) range(rmse, na.rm = T)
# 1            mouse                -0.19                   0.17
# 2            mouse                 0.95                   0.96
# 3         chipmunk                -0.33                   0.00
# 4         chipmunk                 1.00                   0.89
# 5 flying\nsquirrel                -0.48                   0.27
# 6 flying\nsquirrel                 0.66                   0.79

#same as function below
# random %>% dplyr::group_by(sp) %>%
#   summarize(range=range(dhat.random, na.rm=T))

##############
#number of models that converged

converge <- merge(mod.coef.all, data.frame(t=rep(rep(c(0,15,60,1440), times=3),times=4), sp=rep(rep(c("mouse","chipmunk","flying\nsquirrel"), each=4),4), mod=rep(c("SC no.info sep","SC info sep","SC no.info pool","SC info pool"), each=12)), all=T)
converge[is.na(converge$n),]$n <- 0

ggplot(data=converge,
       aes(x=mod, y=n, fill=t)) +geom_hline(aes(yintercept=8), col="grey50",lty="dashed") + geom_bar(stat='identity', position="dodge") + 
  # scale_fill_manual(values=c("#440154FF","#21908CFF","#FDE725FF"), name="model") +
  # scale_fill_manual(values=sc.vals, name="model") +
  ylab("number sites with\nmodel convergence") + xlab("model") +
  scale_y_continuous(breaks=0:8) +
  theme_bw(base_size=16) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5), legend.position = "top") + facet_wrap(~sp)
ggsave(filename="Figures/CAM-SC_fig_nconverged.tiff", height=6, width=8, dpi=400, compression="lzw")



################
#random density ranges for each species
random %>% dplyr::group_by(sp) %>%
  reframe(range=range(dhat.random)) %>% as.data.frame() %>% mutate_if(is.numeric, round, digits=2)

# sp                 range
# <chr>              <dbl>
# 5 "mouse"            13.12 
# 6 "mouse"            62.84 

# 1 "chipmunk"          1.68
# 2 "chipmunk"          8.59

# 3 "flying\nsquirrel"  1.14
# 4 "flying\nsquirrel"  3.02

################
#number of SC models that had positive slopes

mod.coef.scaled[!is.na(mod.coef.scaled$Estimate) & mod.coef.scaled$Estimate > 0,] %>% dplyr::group_by(sp) %>%
  reframe(count=length(Estimate))

# 1 "mouse"                7
# 2 "chipmunk"             7
# 3 "flying\nsquirrel"     3