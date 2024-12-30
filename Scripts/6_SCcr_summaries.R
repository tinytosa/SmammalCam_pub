#########################
#########################

#SCR random on the y axis
#SC from live capture recapture estimates on x axis

#########################
#load packages
require(dplyr)
require(tidyr)

require(ggplot2)
require(ggpubr)
require(viridis)

#########################
#conversion tables

sc.vals <- c("#AAC0AA","#386150","#896A67","#1F1300")

order.mods <- c("CR-SC_noinfo_separate","CR-SC_info_separate","CR-SC_noinfo_pool","CR-SC_info_pool")
order.mods.names <- c("SC no.info sep","SC info sep","SC no.info pool","SC info pool")

order.sp <- c("PEMA","TATO","GLSA")
order.sp.names <- c("mouse","chipmunk","flyingsquirrel","flying\nsquirrel")

#
sp.convert <- data.frame(sp=c("chipmunk","flying\nsquirrel","flying\nsquirrel","mouse"), acro=c("TATO","GLSA","GLOR","PEMA"))
# sp.convert <- data.frame(sp=c("chipmunk","flying\nsquirrel","mouse"), sp2=c("TOWNSENDS CHIPMUNK","FLYING SQUIRREL","DEER MOUSE"))
sp.convert$sp <- factor(sp.convert$sp, levels=c("mouse","chipmunk","flying\nsquirrel"))

#########################
#load data

#random SCR estimates
# random <- read.table(file="Data_linearregs/CR_estimates.txt", sep=",", header=T) #wrong grid numbers for PEMA
random <- read.table(file="Data_linearregs/CR_estimates_correctgridnums.txt", sep=",", header=T)
random <- random[random$model == "SCR random",] #has absolute values and scaled values
random <- random[,c("sp","grid",
                    "dhat","mean.random","sd.random","X2.5..random","X97.5..random", #remove extra columns
                    "dhat.scaled","X2.5.r.scaled","X97.5.r.scaled")]
names(random) <- c("sp","grid", "dhat.random","mean.random","sd.random","X2.5..random","X97.5..random","dhat.random.scaled","X2.5.r.scaled","X97.5.r.scaled")
random[random$sp == "flyingsquirrel",]$sp <- "flying\nsquirrel"

table(random$sp)

sc.cr <- read.table("Output_CR/CR-SC_estimates.txt", sep=",", header=T)
sc.cr <- sc.cr[grep(sc.cr$param, pattern="D"),]
sc.cr$converged <- "no convergence"
sc.cr[sc.cr$Rhat < 1.3,]$converged <- "converged"

table(sc.cr[,c("sp","model")])
table(sc.cr[sc.cr$converged == "converged",c("sp","model")])
# model
# sp   CR-SC_info_pool CR-SC_info_separate CR-SC_noinfo_pool CR-SC_noinfo_separate
# GLSA               7                   8                 0                     5
# PEMA               5                   5                 1                     1
# TATO               9                   8                 9                     4

########################
sc.cr <- sc.cr[sc.cr$converged == "converged",]

sc.cr <- merge(sp.convert, sc.cr, by.y="sp", by.x="acro", all.y=T)
sc.cr <- merge(sc.cr, random, by.x=c("g","sp"), by.y=c("grid","sp"), all.x=T)

#scale sc.cr mean values
#only scale values that converged
sc.cr$mode.scaled <- NA
sc.cr$X2.5.scaled <- NA
sc.cr$X97.5.scaled <- NA

mod.coef.all <- data.frame()
mod.coef.scaled <- data.frame()

  for(sp in unique(sc.cr$sp))
  {
    print(sp)
    for(m in unique(sc.cr$model))
    {
      print(m)
      if(nrow(sc.cr[sc.cr$sp == sp & sc.cr$model == m,]) > 0)
      {
        ########
        #not scaled
        #linear regression with non-scaled values
        l <- lm(dhat.random ~ mode, data=sc.cr[sc.cr$sp == sp & sc.cr$model == m,])
        
        mod.coef <- data.frame(coef(summary(l)))[2,]
        mod.coef$sp <- sp
        mod.coef$mod <- m
        mod.coef$r2 <- summary(l)$adj.r.squared
        mod.coef$rmse <- sqrt(mean(l$residuals^2))
        mod.coef$n <- nrow(sc.cr[sc.cr$sp == sp & sc.cr$model == m,])
        
        mod.coef.all <- rbind(mod.coef.all, mod.coef)
        
        ########
        #scaled
        s <- scale(sc.cr[sc.cr$sp == sp  & sc.cr$model == m,]$mode)
        sc.cr[sc.cr$sp == sp & sc.cr$model == m,]$mode.scaled <- s
        
        sc.cr[sc.cr$sp == sp & sc.cr$model == m,]$X2.5.scaled <- (sc.cr[sc.cr$sp == sp & sc.cr$model == m,]$X2.5. - attr(s, "scaled:center"))/attr(s, "scaled:scale")
        sc.cr[sc.cr$sp == sp & sc.cr$model == m,]$X97.5.scaled <- (sc.cr[sc.cr$sp == sp & sc.cr$model == m,]$X97.5. - attr(s, "scaled:center"))/attr(s, "scaled:scale")
        
        if(nrow(sc.cr[sc.cr$sp == sp & sc.cr$model == m,]) == 1)
        {
          sc.cr[sc.cr$sp == sp & sc.cr$model == m,]$mode.scaled <- 0
        }
        
        l <- lm(dhat.random.scaled ~ mode.scaled, data=sc.cr[sc.cr$sp == sp & sc.cr$model == m,]) #linear regression
        
        mod.coef.s <- data.frame(coef(summary(l)))[2,]
        mod.coef.s$sp <- sp
        mod.coef.s$mod <- m
        mod.coef.s$r2 <- summary(l)$adj.r.squared
        mod.coef.s$rmse <- sqrt(mean(l$residuals^2))
        
        mod.coef.s$n <- nrow(sc.cr[sc.cr$sp == sp & sc.cr$model == m,])
        
        mod.coef.scaled <- rbind(mod.coef.scaled, mod.coef.s)
        
      }
    }
  }

# write.table(sc.cr, file="Data_linearregs/sc-cr_estimates.txt", sep=",", row.names=F)

sc.cr %>% dplyr::group_by(sp, model) %>% reframe(fold = range(mode/dhat.random, na.rm = T)) %>% as.data.frame() %>% mutate_if(is.numeric, round, digits=2)

#                  sp                 model  fold
# 1             mouse       CR-SC_info_pool  3.24
# 2             mouse       CR-SC_info_pool  9.31
# 3             mouse   CR-SC_info_separate 10.17
# 4             mouse   CR-SC_info_separate 23.62
# 5             mouse     CR-SC_noinfo_pool  3.46
# 6             mouse     CR-SC_noinfo_pool  3.46
# 7             mouse CR-SC_noinfo_separate  1.05
# 8             mouse CR-SC_noinfo_separate  1.05

# 9          chipmunk       CR-SC_info_pool  2.87
# 10         chipmunk       CR-SC_info_pool  5.57
# 11         chipmunk   CR-SC_info_separate  6.93
# 12         chipmunk   CR-SC_info_separate 26.50
# 13         chipmunk     CR-SC_noinfo_pool  0.50
# 14         chipmunk     CR-SC_noinfo_pool  0.94
# 15         chipmunk CR-SC_noinfo_separate  0.03
# 16         chipmunk CR-SC_noinfo_separate  0.41

# 17 flying\nsquirrel       CR-SC_info_pool  8.64
# 18 flying\nsquirrel       CR-SC_info_pool 11.67
# 19 flying\nsquirrel   CR-SC_info_separate 19.88
# 20 flying\nsquirrel   CR-SC_info_separate 48.84
# 21 flying\nsquirrel CR-SC_noinfo_separate  0.02
# 22 flying\nsquirrel CR-SC_noinfo_separate  0.29

sc.cr[,c("sp","acro","model","mode","dhat.random")]

mod.coef.all %>% dplyr::group_by(sp) %>%
  reframe(range.slope=range(Estimate, na.rm=T), range.r2=range(r2, na.rm=T), range.rmse=range(rmse, na.rm=T)) %>% as.data.frame() %>% mutate_if(is.numeric, round, digits=2)

#                 sp range.slope range.r2 range.rmse
# 5            mouse        0.11     0.00       0.00
# 6            mouse        0.13     0.42      15.23

# 1         chipmunk       -0.45    -0.49       0.76
# 2         chipmunk        1.52     0.89       2.69

# 3 flying\nsquirrel        0.05     0.17       0.21
# 4 flying\nsquirrel        1.42     0.86       0.51

mod.coef.scaled %>% dplyr::group_by(sp) %>%
  reframe(range.slope=range(Estimate, na.rm=T), range.r2=range(r2, na.rm=T), range.rmse=range(rmse, na.rm=T)) %>% as.data.frame() %>% mutate_if(is.numeric, round, digits=2)
#                 sp range.slope range.r2 range.rmse
# 5            mouse        0.77     0.00       0.00
# 6            mouse        0.79     0.42       0.89

# 1         chipmunk       -0.09    -0.49       0.29
# 2         chipmunk        0.95     0.89       1.04

# 3 flying\nsquirrel        0.50     0.17       0.29
# 4 flying\nsquirrel        0.87     0.86       0.74

mod.coef.scaled[!mod.coef.scaled$mod %in% c("CR-SC_noinfo_separate"),] %>% reframe(range.slope=range(Estimate, na.rm=T), range.r2=range(r2, na.rm=T), range.rmse=range(rmse, na.rm=T)) %>% as.data.frame() %>% mutate_if(is.numeric, round, digits=2)
#   range.slope range.r2 range.rmse
# 1        0.29    -0.06       0.00
# 2        0.95     0.89       0.89

mod.coef.scaled[mod.coef.scaled$mod %in% c("CR-SC_noinfo_separate") & mod.coef.scaled$sp == "flying\nsquirrel",] %>% mutate_if(is.numeric, round, digits=2)
#              Estimate Std..Error t.value Pr...t..               sp                   mod   r2 rmse n
# mode.scaled5     0.73       0.48    1.54     0.22 flying\nsquirrel CR-SC_noinfo_separate 0.26 0.74 5

mod.coef.scaled[mod.coef.scaled$mod %in% c("CR-SC_noinfo_separate") & mod.coef.scaled$sp == "mouse",] %>% mutate_if(is.numeric, round, digits=2)
mod.coef.scaled[mod.coef.scaled$mod %in% c("CR-SC_noinfo_separate") & mod.coef.scaled$sp == "chipmunk",] %>% mutate_if(is.numeric, round, digits=2)
#################
#1b. plot unscaled values together

mod.coef.all$mod <- factor(mod.coef.all$mod, levels=order.mods, labels=order.mods.names)
mod.coef.all$sp <- factor(mod.coef.all$sp, levels=c("mouse","chipmunk","flyingsquirrel","flying\nsquirrel"))

sc.cr$model <- factor(sc.cr$model, levels=order.mods, labels=order.mods.names)

all.0 <- ggplot(data=sc.cr, aes(x=mode, y=dhat.random, col=model, shape=model)) +
  geom_abline(aes(slope=1, intercept=0), col="grey50", lwd=0.5) +
  geom_smooth(method="lm", aes(fill=model, lty=model), alpha=0.1) + #add ribbons
  geom_errorbar(aes(ymin=X2.5..random, ymax=X97.5..random), width=0.01, lty="dashed", linewidth=0.1) + #vertical error bars
  geom_errorbar(aes(xmin=X2.5., xmax=X97.5.), width=0.01, lty="dashed", linewidth=0.1) +
  geom_point(size=3, stroke=2) + #, alpha=0.5) +
  scale_color_manual(values=sc.vals) + scale_fill_manual(values=sc.vals) +
  scale_linetype(name="model") +
  scale_shape_manual(name="model", values=c(2, 2, 1, 1)) +
  xlab("CR SC density estimate (animals/ha)") + ylab("SCR random (animals/ha)") +
  theme_bw(base_size=20) + theme(panel.grid = element_blank(), legend.position="top", strip.background = element_blank(), strip.text.x = element_blank()) + 
  facet_wrap(~sp, scales="free") +
  coord_cartesian(ylim=c(0, NA))
# coord_cartesian(ylim=c(0, 150))
all.0

#slope
s <- ggplot(mod.coef.all, aes(x=mod, y=Estimate, col=mod, shape=mod)) + 
  geom_point(size=5, position=position_dodge(width=1)) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=1)) +
  scale_color_manual(values=sc.vals) +
  scale_shape_manual(name="model", values=c(17, 17, 16, 16)) +
  geom_hline(yintercept=0, lty="dashed") +
  geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("slope") + #scale_y_continuous(limits=c(0,NA)) +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + 
  facet_wrap(~sp, scale="free_y")

#r2
r2 <- ggplot(mod.coef.all, aes(x=mod, y=r2, col=mod, shape=mod)) + geom_point(size=5, position=position_dodge(width=1)) + 
  scale_color_manual(values=sc.vals) +
  scale_shape_manual(name="model", values=c(17, 17, 16, 16)) +
  geom_hline(yintercept = 0, lty="dashed") + geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("R2") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + 
  facet_wrap(~sp, scale="free_y")

#rmse
rmse <- ggplot(mod.coef.all, aes(x=mod, y=rmse, col=mod, shape=mod)) + geom_point(size=5, position=position_dodge(width=1)) + 
  scale_color_manual(values=sc.vals) +
  scale_shape_manual(name="model", values=c(17, 17, 16, 16)) +
  geom_hline(yintercept = 0, lty="dashed") + #geom_hline(yintercept=1, lty="dashed") + 
  xlab("model") + ylab("RMSE") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5)) + facet_wrap(~sp, scale="free_y")

p <- ggarrange(all.0, s, r2, rmse, ncol=1, align="v", heights=c(0.5,0.25,0.15,0.25))
p
# ggsave(p, filename="Figures/CR-SC_fig.tiff", height=15, width=15, units="in", dpi=400, compression="lzw")

################
#1c. scale values and plot and calculate slopes, RMSE, and R2
mod.coef.scaled$sp <- factor(mod.coef.scaled$sp, levels=c("mouse","chipmunk","flyingsquirrel","flying\nsquirrel"))
mod.coef.scaled$mod <- factor(mod.coef.scaled$mod, levels=order.mods, labels=order.mods.names)

all.0 <- ggplot(data=sc.cr, aes(x=mode.scaled, y=dhat.random, col=model, shape=model)) +
  geom_smooth(data=sc.cr, method="lm", aes(fill=model, lty=model), alpha=0.1) + #add ribbons
  geom_errorbar(aes(ymin=X2.5..random, ymax=X97.5..random), width=0.01, lty="dashed", linewidth=0.1) + #vertical error bars
  geom_errorbar(aes(xmin=X2.5.scaled, xmax=X97.5.scaled), width=0.01, lty="dashed", linewidth=0.1) +
  geom_point(size=3, stroke=1) + #, alpha=0.5) +
  scale_color_manual(values=sc.vals) + scale_fill_manual(values=sc.vals) +
  # scale_shape_manual(name="model", values=c(17, 17, 16, 16)) +
  scale_linetype(name="model") +
  scale_shape_manual(name="model", values=c(2, 2, 1, 1)) +
  xlab("scaled model estimate (detections/camera)") + ylab("SCR random (animals/ha)") +
  theme_bw(base_size=20) + theme(panel.grid = element_blank(), legend.position="top", strip.background = element_blank(), strip.text.x = element_blank()) + 
  facet_wrap(~sp, scales="free") +
  coord_cartesian(ylim=c(0, NA), xlim=c(NA, 3))
# coord_cartesian(ylim=c(0, 150))
all.0

#slope
s <- ggplot(mod.coef.scaled, aes(x=mod, y=Estimate, col=mod, shape=mod)) + 
  geom_point(size=5, position=position_dodge(width=0.75)) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=0.75)) +
  scale_color_manual(values=sc.vals) +
  scale_shape_manual(name="model", values=c(17, 17, 16, 16)) +
  geom_hline(yintercept=0, lty="dashed") +
  geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("slope") + #scale_y_continuous(limits=c(0,NA)) +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + 
  facet_wrap(~sp, scale="free_y")

#r2
r2 <- ggplot(mod.coef.scaled, aes(x=mod, y=r2, col=mod, shape=mod)) + geom_point(size=5, position=position_dodge(width=0.75)) + 
  scale_color_manual(values=sc.vals) +
  scale_shape_manual(name="model", values=c(17, 17, 16, 16)) +
  geom_hline(yintercept = 0, lty="dashed") + geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("R2") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + 
  facet_wrap(~sp, scale="free_y")

#rmse
rmse <- ggplot(mod.coef.scaled, aes(x=mod, y=rmse, col=mod, shape=mod)) + geom_point(size=5, position=position_dodge(width=0.75)) + 
  scale_color_manual(values=sc.vals) +
  scale_shape_manual(name="model", values=c(17, 17, 16, 16)) +
  geom_hline(yintercept = 0, lty="dashed") + #geom_hline(yintercept=1, lty="dashed") + 
  xlab("model") + ylab("RMSE") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5)) + facet_wrap(~sp, scale="free_y")

p <- ggarrange(all.0, s, r2, rmse, ncol=1, align="v", heights=c(0.5,0.25,0.15,0.25))
p
# ggsave(p, filename="Figures/CR-SC_scaled_fig.tiff", height=15, width=15, units="in", dpi=400, compression="lzw")

write.csv(mod.coef.scaled, file="Data_linearregs/CR-SC_slopes.csv")

detach("package:plyr", unload = TRUE)
mod.coef.scaled %>% group_by(sp) %>%
  summarize(min(Estimate, na.rm=T), max(Estimate, na.rm=T), max(r2), min(rmse))

#   sp                 `min(Estimate, na.rm = T)` `max(Estimate, na.rm = T)` `max(r2)` `min(rmse)`
#   <fct>                                   <dbl>                      <dbl>     <dbl>       <dbl>
# 1 "mouse"                                0.775                       0.792     0.422       0    
# 2 "chipmunk"                           -0.0910                       0.950     0.888       0.295
# 3 "flying\nsquirrel"                     0.498                       0.874     0.860       0.295


#############
#converged?
ggplot(data=rbind(mod.coef.all, data.frame(Estimate=0, Std..Error=NA, t.value=NA, Pr...t..=NA, sp="flying\nsquirrel", mod="SC no.info pool", r2=NA, rmse=NA, n=0)),
       aes(x=mod, y=n, fill=mod)) +geom_hline(aes(yintercept=9), col="grey50",lty="dashed") + geom_bar(stat='identity', position="dodge") + 
  # scale_fill_manual(values=c("#440154FF","#21908CFF","#FDE725FF"), name="model") +
  scale_fill_manual(values=sc.vals, name="model") +
  ylab("number sites with\nmodel convergence") + xlab("model") +
  scale_y_continuous(breaks=0:9) +
  theme_bw(base_size=16) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5), legend.position = "none") + facet_wrap(~sp)
ggsave(filename="Figures/CR-SC_fig_nconverged.tiff", height=4, width=8, dpi=400, compression="lzw")
