#########################
#########################
#data represents number of detections per camera in sherman grid for PEMA and tomahawk grid in TATO and GLOR
#SCR random on the y axis
#Nmixture estimates on x axis

#convert to density

#########################
#load packages
require(dplyr)

require(ggplot2)
require(ggpubr)
require(viridis)

#########################
#load data

#conversion table
sp.convert <- data.frame(sp=c("chipmunk","flying\nsquirrel","flying\nsquirrel","mouse"), acro=c("TATO","GLSA","GLOR","PEMA"))
sp.convert$sp <- factor(sp.convert$sp, levels=c("mouse","chipmunk","flying\nsquirrel"))

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

#random SCR estimates
# random <- read.table(file="Data_linearregs/CR_estimates.txt", sep=",", header=T) #wrong grid numbers for PEMA
random <- read.table(file="Data_linearregs/CR_estimates_correctgridnums.txt", sep=",", header=T)
random <- random[random$model == "SCR random",] #has absolute values and scaled values
random <- random[,c("sp","grid",
                    "dhat","mean.random","sd.random","X2.5..random","X97.5..random", #remove extra columns
                    "dhat.scaled","X2.5.r.scaled","X97.5.r.scaled")]
names(random) <- c("sp","grid", "dhat.random","mean.random","sd.random","X2.5..random","X97.5..random","dhat.random.scaled","X2.5.r.scaled","X97.5.r.scaled")
table(random$sp)

#Nmixture model estimates
nmix <- read.table("Output_cam//output_Nmix_2023-09-26.txt", sep=",")
nmix <- nmix[nmix$param == "lambda",]
nmix <- nmix[nmix$grid > 0,]
nmix <- nmix[nmix$sp != "myo",]

nmix$sp <- toupper(nmix$sp)

#turn values into factors
nmix <- merge(sp.convert, nmix, by.y="sp", by.x="acro", all.y=T) #factor species names
nmix$t <- factor(nmix$ctime, levels=c(0,15,60,1440)) #factor delta 
nmix$mod <- factor(nmix$mod, levels=c("base","psite","pdecay"), labels=c("base","pstation","pdecay")) #factor model

nmix <- merge(nmix, random, by=c("grid","sp"))

nmix <- merge(nmix, info[,c("sp","mmdm.a","halfmmdm.a")], by.x="acro", by.y="sp", all.x=T)
nmix$mean.density <- nmix$mean/nmix$halfmmdm.a

#scale nmix mean values
#only scale values that converged
nmix$mean.scaled <- NA
nmix$X2.5.scaled <- NA
nmix$X97.5.scaled <- NA

mod.coef.all <- data.frame()
mod.coef.scaled <- data.frame()

for(t in unique(nmix$ctime))
{
  print(t)
  for(sp in unique(nmix$sp))
  {
    print(sp)
    for(m in unique(nmix$mod))
    {
      print(m)
      if(nrow(nmix[nmix$ctime == t & nmix$sp == sp & nmix$converged == "converged" & nmix$mod == m,]) > 0)
      {
        ########
        #not scaled
        #linear regression with non-scaled values
        # l <- lm(dhat.random ~ mean, data=nmix[nmix$ctime == t & nmix$sp == sp & nmix$converged == "converged" & nmix$mod == m,])
        l <- lm(dhat.random ~ mean/halfmmdm.a, data=nmix[nmix$ctime == t & nmix$sp == sp & nmix$converged == "converged" & nmix$mod == m,])
        
        mod.coef <- data.frame(coef(summary(l)))[2,]
        mod.coef$sp <- sp
        mod.coef$t <- t
        mod.coef$mod <- m
        mod.coef$r2 <- summary(l)$adj.r.squared
        mod.coef$rmse <- sqrt(mean(l$residuals^2))
        mod.coef$n <- nrow(nmix[nmix$ctime == t & nmix$sp == sp & nmix$converged == "converged" & nmix$mod == m,])
        mod.coef.all <- rbind(mod.coef.all, mod.coef)
        
        ########
        #scaled
        s <- scale(nmix[nmix$ctime == t & nmix$sp == sp  & nmix$converged == "converged" & nmix$mod == m,]$mean)
        nmix[nmix$ctime == t & nmix$sp == sp & nmix$converged == "converged" & nmix$mod == m,]$mean.scaled <- s
        
        nmix[nmix$ctime == t & nmix$sp == sp & nmix$converged == "converged" & nmix$mod == m,]$X2.5.scaled <- (nmix[nmix$ctime == t & nmix$sp == sp & nmix$converged == "converged" & nmix$mod == m,]$X2.5. - attr(s, "scaled:center"))/attr(s, "scaled:scale")
        nmix[nmix$ctime == t & nmix$sp == sp & nmix$converged == "converged" & nmix$mod == m,]$X97.5.scaled <- (nmix[nmix$ctime == t & nmix$sp == sp & nmix$converged == "converged" & nmix$mod == m,]$X97.5. - attr(s, "scaled:center"))/attr(s, "scaled:scale")
        
        if(nrow(nmix[nmix$ctime == t & nmix$sp == sp & nmix$converged == "converged" & nmix$mod == m,]) == 1)
        {
          nmix[nmix$ctime == t & nmix$sp == sp & nmix$converged == "converged" & nmix$mod == m,]$mean.scaled <- 0
        }
        
        l <- lm(dhat.random.scaled ~ mean.scaled, data=nmix[nmix$ctime == t & nmix$sp == sp & nmix$converged == "converged" & nmix$mod == m,]) #linear regression
        
        mod.coef.s <- data.frame(coef(summary(l)))[2,]
        mod.coef.s$sp <- sp
        mod.coef.s$t <- t
        mod.coef.s$mod <- m
        mod.coef.s$r2 <- summary(l)$adj.r.squared
        mod.coef.s$rmse <- sqrt(mean(l$residuals^2))
        mod.coef.s$n <- nrow(nmix[nmix$ctime == t & nmix$sp == sp & nmix$converged == "converged" & nmix$mod == m,])
        mod.coef.scaled <- rbind(mod.coef.scaled, mod.coef.s)
      }
    }
  }
}

# write.table(nmix, file="Data_linearregs/nmix_estimates.txt", sep=",", row.names=F)


#####
######
#1b. plot unscaled values together
# t <- 15
# sp <- "mouse"

mod.coef.all$sp <- factor(mod.coef.all$sp, levels=c("mouse","chipmunk","flyingsquirrel","flying\nsquirrel"))
mod.coef.all$t <- factor(mod.coef.all$t, levels=c(0,15,60,1440))
mod.coef.all$mod <- factor(mod.coef.all$mod, levels=c("base","pstation","pdecay"))

# for(t in c(0, 15, 60, 1440)) #t in minutes
# {
#   print(t)
# avgdet.t <- avgdets.scaled[avgdets.scaled$sp == sp & avgdets.scaled$t == t,]
# avgdet.t <- avgdets.scaled[avgdets.scaled$t == t,]
# mod.coef.t <- mod.coef.all[mod.coef.all$t==t,]

# all.0 <- ggplot(data=nmix[nmix$t == 0 & nmix$converged == "converged",], aes(x=mean, y=dhat.random, col=mod, group=mod)) +
all.0 <- ggplot(data=nmix[nmix$t == 0 & nmix$converged == "converged",], aes(x=mean/halfmmdm.a, y=dhat.random, col=mod, group=mod)) +
  geom_smooth(method="lm", aes(fill=mod, lty=mod), alpha=0.1) + #add ribbons
  geom_errorbar(aes(ymin=X2.5..random, ymax=X97.5..random), width=0.01, lty="dashed", linewidth=0.1) + #vertical error bars
  geom_errorbar(aes(xmin=(mean-sd)/halfmmdm.a, xmax=(mean+sd)/halfmmdm.a), width=0.01, lty="dashed", linewidth=0.1) +
  geom_point(size=3, pch=1, stroke=2) + #, alpha=0.5) +
  # scale_color_manual(values=c.vals) + scale_fill_manual(values=c.vals) +
  # scale_color_brewer(palette="BrBG", name="days active") +
  # scale_fill_brewer(palette="BrBG", name="days active") +
  scale_color_viridis(discrete=T, name="model") +
  scale_fill_viridis(discrete=T, name="model") +
  scale_linetype(name="model") +
  xlab("N-mixture density estimate (animals/ha)") + ylab("SCR random (animals/ha)") +
  theme_bw(base_size=20) + theme(panel.grid = element_blank(), legend.position="top", strip.background = element_blank(), strip.text.x = element_blank()) + 
  facet_wrap(~sp, scales="free") +
  coord_cartesian(ylim=c(0, NA))
  # coord_cartesian(ylim=c(0, 150))
all.0

#slope
s <- ggplot(mod.coef.all, aes(x=t, y=Estimate, col=mod, shape=t)) +
# s.main <- ggplot(mod.coef.all[mod.coef.all$t != 1440,], aes(x=t, y=Estimate, col=mod, shape=t)) +
  geom_point(size=5, position=position_dodge(width=1)) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=1)) +
  scale_color_viridis(discrete=T, name="model") +
  geom_hline(yintercept=0, lty="dashed") +
  geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("slope") + #scale_y_continuous(limits=c(0,NA)) +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + 
  facet_wrap(~sp, scale="free_y")

#r2
r2 <- ggplot(mod.coef.all, aes(x=t, y=r2, col=mod, shape=t)) + geom_point(size=5, position=position_dodge(width=1)) + 
# r2 <- ggplot(mod.coef.all[mod.coef.all$t != 1440,], aes(x=t, y=r2, col=mod, shape=t)) + geom_point(size=5, position=position_dodge(width=1)) + 
  # scale_color_manual(values=c.vals) +
  scale_color_viridis(discrete=T, name="model") +
  geom_hline(yintercept = 0, lty="dashed") + geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("R2") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + 
  facet_wrap(~sp, scale="free_y")

#rmse
rmse <- ggplot(mod.coef.all, aes(x=t, y=rmse, col=mod, shape=t)) + geom_point(size=5, position=position_dodge(width=1)) + 
# rmse <- ggplot(mod.coef.all[mod.coef.all$t != 1440,], aes(x=t, y=rmse, col=mod, shape=t)) + geom_point(size=5, position=position_dodge(width=1)) + 
  # scale_color_manual(values=c.vals) +
  scale_color_viridis(discrete=T, name="model") +
  geom_hline(yintercept = 0, lty="dashed") + #geom_hline(yintercept=1, lty="dashed") + 
  xlab("model") + ylab("RMSE") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5)) + facet_wrap(~sp, scale="free_y")

p <- ggarrange(all.0, s, r2, rmse, ncol=1, align="v", heights=c(0.5,0.25,0.15,0.25))
p
# ggsave(p, filename="Figures/nmix_fig.tiff", height=15, width=15, units="in", dpi=400, compression="lzw")
# ggsave(p, filename="Figures/nmix_fig_no1440.tiff", height=15, width=15, units="in", dpi=400, compression="lzw")

# }

#########################
#1c. scale values and plot and calculate slopes, RMSE, and R2
mod.coef.scaled$sp <- factor(mod.coef.scaled$sp, levels=c("mouse","chipmunk","flyingsquirrel","flying\nsquirrel"))
mod.coef.scaled$t <- factor(mod.coef.scaled$t, levels=c(0,15,60,1440))
mod.coef.scaled$mod <- factor(mod.coef.scaled$mod, levels=c("base","pstation","pdecay"))

all.0 <- ggplot(data=nmix[nmix$t == 0,], aes(x=mean.scaled, y=dhat.random, col=mod, group=mod)) +
  geom_smooth(data=nmix[nmix$t == 0,], method="lm", aes(fill=mod, lty=mod), alpha=0.1) + #add ribbons
  geom_errorbar(aes(ymin=X2.5..random, ymax=X97.5..random), width=0.01, lty="dashed", linewidth=0.1) + #vertical error bars
  geom_errorbar(aes(xmin=X2.5.scaled, xmax=X97.5.scaled), width=0.01, lty="dashed", linewidth=0.1) +
  geom_point(size=3, pch=1, stroke=2) + #, alpha=0.5) +
  scale_color_viridis(discrete=T, name="model") +
  scale_fill_viridis(discrete=T, name="model") +
  scale_linetype(name="model") +
  xlab("scaled model estimate (detections/camera)") + ylab("SCR random (animals/ha)") +
  theme_bw(base_size=20) + theme(panel.grid = element_blank(), legend.position="top", strip.background = element_blank(), strip.text.x = element_blank()) + 
  facet_wrap(~sp, scales="free") +
  coord_cartesian(ylim=c(0, NA))
  # coord_cartesian(ylim=c(0, 150))
all.0

#slope
s <- ggplot(mod.coef.scaled, aes(x=t, y=Estimate, col=mod, shape=t)) +
# s <- ggplot(mod.coef.scaled[mod.coef.scaled$t != 1440,], aes(x=t, y=Estimate, col=mod, shape=t)) + 
  geom_point(size=5, position=position_dodge(width=0.75)) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=0.75)) +
  scale_color_viridis(discrete=T, name="model") +
  geom_hline(yintercept=0, lty="dashed") +
  geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("slope") + #scale_y_continuous(limits=c(0,NA)) +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + 
  facet_wrap(~sp, scale="free_y")

#r2
r2 <- ggplot(mod.coef.scaled[mod.coef.scaled$n > 1,], aes(x=t, y=r2, col=mod, shape=t)) + geom_point(size=5, position=position_dodge(width=0.75)) +
# r2 <- ggplot(mod.coef.scaled[mod.coef.scaled$n > 1 & mod.coef.scaled$t != 1440,], aes(x=t, y=r2, col=mod, shape=t)) + geom_point(size=5, position=position_dodge(width=0.75)) + 
  scale_color_viridis(discrete=T, name="model") +
  geom_hline(yintercept = 0, lty="dashed") + geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("R2") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + 
  facet_wrap(~sp, scale="free_y")

#rmse
rmse <- ggplot(mod.coef.scaled[mod.coef.scaled$n > 1,], aes(x=t, y=rmse, col=mod, shape=t)) + geom_point(size=5, position=position_dodge(width=0.75)) +
# rmse <- ggplot(mod.coef.scaled[mod.coef.scaled$n > 1 & mod.coef.scaled$t != 1440,], aes(x=t, y=rmse, col=mod, shape=t)) + geom_point(size=5, position=position_dodge(width=0.75)) + 
  scale_color_viridis(discrete=T, name="model") +
  geom_hline(yintercept = 0, lty="dashed") + #geom_hline(yintercept=1, lty="dashed") + 
  xlab("model") + ylab("RMSE") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5)) + facet_wrap(~sp, scale="free_y")

p <- ggarrange(all.0, s, r2, rmse, ncol=1, align="v", heights=c(0.5,0.25,0.15,0.25))
p
# ggsave(p, filename="Figures/nmix_scaled_fig.tiff", height=15, width=15, units="in", dpi=400, compression="lzw")
# ggsave(p, filename="Figures/nmix_scaled_fig_2.tiff", height=15, width=15, units="in", dpi=400, compression="lzw") #to show mouse top panel

write.csv(mod.coef.scaled, file="Data_linearregs/nmix_slopes.csv")

############
detach("package:plyr", unload = TRUE)
mod.coef.scaled %>% group_by(sp) %>%
  summarize(min(Estimate, na.rm=T), max(Estimate, na.rm=T), max(r2), min(rmse)) %>% mutate_if(is.numeric, round, digits=2)
#    sp                 `min(Estimate, na.rm = T)` `max(Estimate, na.rm = T)` `max(r2)` `min(rmse)`
# <fct>                                   <dbl>                      <dbl>     <dbl>      <dbl>
# 1 "mouse"                               -0.303                       0.86   NaN         0     
# 2 "chipmunk"                            -0.614                       1.18     0.99      0.05
# 3 "flying\nsquirrel"                    -0.013                       0.59     0.48      0.42 

mod.coef.scaled[order(mod.coef.scaled$sp, mod.coef.scaled$mod, mod.coef.scaled$t),] %>% mutate_if(is.numeric, round, 3)
#               Estimate Std..Error t.value Pr...t..               sp    t      mod     r2  rmse n
# mean.scaled14    0.181      0.409   0.443    0.681            mouse    0     base -0.192 0.748 6
# mean.scaled31    0.042      0.525   0.081    0.941            mouse   15     base -0.330 0.813 5
# mean.scaled5    -0.051      0.463  -0.110    0.917            mouse   60     base -0.197 0.958 7
# mean.scaled22    0.169      0.457   0.371    0.726            mouse 1440     base -0.168 0.946 7
# mean.scaled12    0.689        NaN     NaN      NaN            mouse    0 pstation    NaN 0.000 2
# mean.scaled29    0.864      0.255   3.384    0.020            mouse   15 pstation  0.635 0.529 7
# mean.scaled3    -0.303      0.443  -0.683    0.525            mouse   60 pstation -0.098 0.917 7
# mean.scaled20    0.241      0.451   0.535    0.615            mouse 1440 pstation -0.135 0.933 7
# mean.scaled13   -0.181        NaN     NaN      NaN            mouse    0   pdecay    NaN 0.000 2
# mean.scaled30    0.498      0.406   1.225    0.275            mouse   15   pdecay  0.077 0.841 7
# mean.scaled4    -0.127      0.460  -0.277    0.793            mouse   60   pdecay -0.182 0.952 7
# mean.scaled21    0.217      0.453   0.479    0.652            mouse 1440   pdecay -0.147 0.938 7

# mean.scaled34    0.925      0.057  16.224    0.039         chipmunk   15     base  0.992 0.047 3
# mean.scaled8     0.427      0.448   0.952    0.442         chipmunk   60     base -0.032 0.549 4
# mean.scaled25   -0.614      0.275  -2.232    0.067         chipmunk 1440     base  0.363 0.631 8
# mean.scaled15    1.180      0.201   5.869    0.028         chipmunk    0 pstation  0.918 0.246 4
# mean.scaled32    0.613      0.276   2.225    0.068         chipmunk   15 pstation  0.361 0.632 8
# mean.scaled6     0.545      0.299   1.823    0.118         chipmunk   60 pstation  0.249 0.684 8
# mean.scaled23   -0.312      0.350  -0.890    0.408         chipmunk 1440 pstation -0.031 0.802 8
# mean.scaled16    0.437      0.327   1.337    0.230         chipmunk    0   pdecay  0.101 0.749 8
# mean.scaled33    0.542      0.300   1.808    0.121         chipmunk   15   pdecay  0.245 0.687 8
# mean.scaled7     0.620      0.273   2.271    0.064         chipmunk   60   pdecay  0.373 0.626 8
# mean.scaled24   -0.519      0.306  -1.694    0.141         chipmunk 1440   pdecay  0.211 0.702 8

# mean.scaled11    0.585      0.271   2.158    0.120 flying\nsquirrel    0     base  0.478 0.420 5
# mean.scaled28    0.064      0.376   0.169    0.871 flying\nsquirrel   15     base -0.161 0.861 8
# mean.scaled2     0.108      0.374   0.288    0.783 flying\nsquirrel   60     base -0.151 0.857 8
# mean.scaled19    0.195      0.368   0.530    0.615 flying\nsquirrel 1440     base -0.115 0.844 8
# mean.scaled9     0.467      0.325   1.437    0.201 flying\nsquirrel    0 pstation  0.132 0.745 8
# mean.scaled26    0.099      0.375   0.264    0.800 flying\nsquirrel   15 pstation -0.153 0.858 8
# mean.scaled      0.112      0.374   0.299    0.775 flying\nsquirrel   60 pstation -0.150 0.857 8
# mean.scaled17    0.221      0.366   0.603    0.568 flying\nsquirrel 1440 pstation -0.100 0.838 8
# mean.scaled10    0.421      0.335   1.254    0.257 flying\nsquirrel    0   pdecay  0.075 0.769 8
# mean.scaled27   -0.013      0.377  -0.034    0.974 flying\nsquirrel   15   pdecay -0.166 0.863 8
# mean.scaled1     0.036      0.377   0.095    0.927 flying\nsquirrel   60   pdecay -0.165 0.863 8
# mean.scaled18    0.164      0.371   0.442    0.674 flying\nsquirrel 1440   pdecay -0.130 0.850 8

##################

#plot each species separately and then combine at the end
main <- c(0, 15, 60, 1440)
sub <- c(0,15,60)
for(sp in unique(mod.coef.all$sp))
{
  print(sp)
  assign(paste("all.", gsub(sp, pattern="\n", replacement=""), sep=""), 
         ggplot(data=nmix[nmix$sp == sp & nmix$t == 0 & nmix$converged == "converged",], aes(x=mean/halfmmdm.a, y=dhat.random, col=mod, group=mod)) +
           geom_smooth(method="lm", aes(fill=mod), alpha=0.1) + #add ribbons #, lty=mod
           geom_errorbar(aes(ymin=X2.5..random, ymax=X97.5..random), width=0.01, lty="dashed", linewidth=0.1) + #vertical error bars
           geom_errorbar(aes(xmin=(mean-sd)/halfmmdm.a, xmax=(mean+sd)/halfmmdm.a), width=0.01, lty="dashed", linewidth=0.1) +
           geom_point(size=3, pch=1, stroke=2) + #, alpha=0.5) +
           # scale_color_manual(values=c.vals) + scale_fill_manual(values=c.vals) +
           # scale_color_brewer(palette="BrBG", name="days active") +
           # scale_fill_brewer(palette="BrBG", name="days active") +
           scale_color_viridis(discrete=T, name="model") +
           scale_fill_viridis(discrete=T, name="model") +
           scale_linetype(name="model") +
           # xlab("N-mixture density estimate (animals/ha)") + ylab("SCR random (animals/ha)") +
           theme_bw(base_size=20) + theme(panel.grid = element_blank(), legend.position="top", strip.background = element_blank(), strip.text.x = element_blank()) + 
           theme(axis.title=element_blank()) + #remove all axis titles
           # facet_wrap(~sp, scales="free") +
           coord_cartesian(ylim=c(0, NA)))
  #slope
  assign(paste("s.", gsub(sp, pattern="\n", replacement=""), ".main", sep=""),
         ggplot(mod.coef.all[mod.coef.all$t %in% main & mod.coef.all$sp == sp,], aes(x=t, y=Estimate, col=mod, shape=t)) + 
           geom_point(size=5, position=position_dodge(width=1)) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=1)) +
           scale_color_viridis(discrete=T, name="model") +
           geom_hline(yintercept=0, lty="dashed") +
           geom_hline(yintercept=1, lty="dashed", col="grey50") + 
           xlab("") + ylab("slope") + #scale_y_continuous(limits=c(0,NA)) +
           theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank(), axis.title=element_blank()) #+
           # facet_wrap(~sp, scale="free_y")
         )
  assign(paste("s.", gsub(sp, pattern="\n", replacement=""), ".sub", sep=""),
         ggplot(mod.coef.all[mod.coef.all$t %in% sub & mod.coef.all$sp == sp,], aes(x=t, y=Estimate, col=mod, shape=t)) + 
           geom_point(size=5, position=position_dodge(width=1)) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=1)) +
           scale_color_viridis(discrete=T, name="model") +
           geom_hline(yintercept=0, lty="dashed") +
           geom_hline(yintercept=1, lty="dashed", col="grey50") + 
           xlab("") + ylab("") + #scale_y_continuous(limits=c(0,NA)) +
           theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), plot.margin = unit(c(0, 0, 0, 0),"inches")) #+ #, axis.text.x = element_blank(), axis.title=element_blank()
           # facet_wrap(~sp, scale="free_y")
         )
  assign(paste("r2.", gsub(sp, pattern="\n", replacement=""), sep=""),
         ggplot(mod.coef.all[mod.coef.all$sp == sp,], aes(x=t, y=r2, col=mod, shape=t)) + geom_point(size=5, position=position_dodge(width=1)) + 
           # scale_color_manual(values=c.vals) +
           scale_color_viridis(discrete=T, name="model") +
           geom_hline(yintercept = 0, lty="dashed") + geom_hline(yintercept=1, lty="dashed") + 
           xlab("") + ylab("R2") +
           theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank(), axis.title=element_blank())) #+ 
  # facet_wrap(~sp, scale="free_y"))
  assign(paste("rmse.", gsub(sp, pattern="\n", replacement=""), sep=""),
         ggplot(mod.coef.all[mod.coef.all$sp == sp,], aes(x=t, y=rmse, col=mod, shape=t)) + geom_point(size=5, position=position_dodge(width=1)) + 
           # scale_color_manual(values=c.vals) +
           scale_color_viridis(discrete=T, name="model") +
           geom_hline(yintercept = 0, lty="dashed") + #geom_hline(yintercept=1, lty="dashed") + 
           xlab("model") + ylab("RMSE") +
           theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5), axis.title=element_blank())) #+ facet_wrap(~sp, scale="free_y"))
}

all <- ggarrange(all.mouse, all.chipmunk, all.flyingsquirrel, nrow=1, common.legend = T)
all <- annotate_figure(all, left=text_grob("SCR random (animals/ha)", rot=90, vjust=1, size=20), bottom=text_grob("N-mixture density estimate (animals/ha)", size=20))

# s.mouse <- ggarrange(s.mouse.t.0, s.mouse.t.15, s.mouse.t.60, s.mouse.t.1440, nrow=1, common.legend = T)
# s.chipmunk <- ggarrange(s.chipmunk.t.0, s.chipmunk.t.15, s.chipmunk.t.60, s.chipmunk.t.1440, nrow=1, common.legend = T)
# s.flyingsquirrel <- ggarrange(s.flyingsquirrel.t.0, s.flyingsquirrel.t.15, s.flyingsquirrel.t.60, s.flyingsquirrel.t.1440, nrow=1, common.legend = T)
# s <- ggarrange(s.mouse.t.0, s.mouse.t.15, s.mouse.t.60, s.mouse.t.1440,
#                s.chipmunk.t.0, s.chipmunk.t.15, s.chipmunk.t.60, s.chipmunk.t.1440,
#                s.flyingsquirrel.t.0, s.flyingsquirrel.t.15, s.flyingsquirrel.t.60, s.flyingsquirrel.t.1440,
#                nrow=1, common.legend = T)
# s <- annotate_figure(s, left=text_grob("slope", rot=90, vjust=1, size=20))
# 
# r2 <- ggarrange(r2.mouse.t.0, r2.mouse.t.15, r2.mouse.t.60, r2.mouse.t.1440,
#                 r2.chipmunk.t.0, r2.chipmunk.t.15, r2.chipmunk.t.60, r2.chipmunk.t.1440,
#                 r2.flyingsquirrel.t.0, r2.flyingsquirrel.t.15, r2.flyingsquirrel.t.60, r2.flyingsquirrel.t.1440,
#                nrow=1, common.legend = T)
# r2 <- annotate_figure(r2, left=text_grob("R2", rot=90, vjust=1, size=20))
# 
# rmse <- ggarrange(rmse.mouse.t.0 + rremove("ylab"), rmse.mouse.t.15 + rremove("ylab"), rmse.mouse.t.60 + rremove("ylab"), rmse.mouse.t.1440 + rremove("ylab"),
#                   rmse.chipmunk.t.0, rmse.chipmunk.t.15, rmse.chipmunk.t.60, rmse.chipmunk.t.1440,
#                   rmse.flyingsquirrel.t.0, rmse.flyingsquirrel.t.15, rmse.flyingsquirrel.t.60, rmse.flyingsquirrel.t.1440,
#                nrow=1, common.legend = T, align="h", widths=rep(0.5,12))
# rmse <- annotate_figure(rmse, left=text_grob("RMSE", rot=90, vjust=1, size=20), bottom=text_grob("consolidation time", size=20))

# main <- c(0,15)
sub <- c(0,15)
# sub <- c(0,15,60)
sp <- "mouse"
# s.mouse.main <- ggplot(mod.coef.all[mod.coef.all$t %in% main & mod.coef.all$sp == sp,], aes(x=t, y=Estimate, col=mod, shape=t)) + 
#   geom_point(size=5, position=position_dodge(width=1)) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=1)) +
#   scale_color_viridis(discrete=T, name="model") +
#   geom_hline(yintercept=0, lty="dashed") +
#   geom_hline(yintercept=1, lty="dashed") + 
#   xlab("") + ylab("slope") + #scale_y_continuous(limits=c(0,NA)) +
#   theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank(), axis.title=element_blank()) + 
#   facet_wrap(~sp, scale="free_y")
s.mouse.sub <- ggplot(mod.coef.all[mod.coef.all$t %in% sub & mod.coef.all$sp == sp,], aes(x=t, y=Estimate, col=mod, shape=t)) + 
  geom_point(size=5, position=position_dodge(width=1)) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=1)) +
  scale_color_viridis(discrete=T, name="model") +
  geom_hline(yintercept=0, lty="dashed") +
  geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("") + #scale_y_continuous(limits=c(0,NA)) +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), plot.margin = unit(c(0, 0, 0, 0),"inches")) #+ #, axis.text.x = element_blank(), axis.title=element_blank()
  # facet_wrap(~sp, scale="free_y")

# main <- c(0,15,60)
# sp <- "chipmunk"
# s.chipmunk.main <- ggplot(mod.coef.all[mod.coef.all$t %in% main & mod.coef.all$sp == sp,], aes(x=t, y=Estimate, col=mod, shape=t)) + 
#   geom_point(size=5, position=position_dodge(width=1)) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=1)) +
#   scale_color_viridis(discrete=T, name="model") +
#   geom_hline(yintercept=0, lty="dashed") +
#   geom_hline(yintercept=1, lty="dashed") + 
#   xlab("") + ylab("slope") + #scale_y_continuous(limits=c(0,NA)) +
#   theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank(), axis.title=element_blank()) + 
#   facet_wrap(~sp, scale="free_y")

require(scales)
show_col(viridis(3))
# s.mouse.sub <- s.mouse.sub + scale_shape_manual(values=c(15,3))
all.chipmunk <- all.chipmunk + scale_color_manual(values=c("#21908CFF","#FDE725FF"), "model") + scale_fill_manual(values=c("#21908CFF","#FDE725FF"), "model") #+ scale_linetype_manual(values=c(2,3), "model")

#combine by row
# s <- ggarrange(s.mouse.main, s.mouse.sub, s.chipmunk.main, s.chipmunk.sub, s.flyingsquirrel.main, nrow=1, widths = c(2,2.25,3,1.5,4), align="h")
# s <- annotate_figure(s, left=text_grob("slope", rot=90, vjust=1, size=20))

s.mouse <- s.mouse.main + annotation_custom(ggplotGrob(s.mouse.sub+rremove("ylab")+rremove("xlab")), xmin=0.5, xmax=2.5, ymin=10, ymax=160)
s.chipmunk <- s.chipmunk.main + annotation_custom(ggplotGrob(s.chipmunk.sub+rremove("ylab")+rremove("xlab")), xmin=0.5, xmax=3.5, ymin=-23, ymax=-5)
s.flyingsquirrel <- s.flyingsquirrel.main + annotation_custom(ggplotGrob(s.flyingsquirrel.sub), xmin=0.5, xmax=3.5, ymin=0.5, ymax=2.25)
# 
# r2 <- ggarrange(r2.mouse, r2.chipmunk, r2.flyingsquirrel, nrow=1, align="h")
# r2 <- annotate_figure(r2, left=text_grob("R2", rot=90, vjust=1, size=20))
# 
# rmse <- ggarrange(rmse.mouse, rmse.chipmunk, rmse.flyingsquirrel, nrow=1, align="h")
# rmse <- annotate_figure(rmse, left=text_grob("RMSE", rot=90, vjust=1, size=20), bottom=text_grob("consolidation time (min)", size=20))
# 
# abs <- ggarrange(all, s, r2, rmse, ncol=1, align="v", heights = c(1,0.75,0.5,0.75))
# 
# ggsave(abs, filename="Figures/nmix_fig_sep.tiff", height=20, width=20, units="in", dpi=400, compression="lzw")

#combine in grid
require(cowplot)
# plot_grid(arrangeGrob(all.mouse, all.chipmunk, all.flyingsquirrel, nrow=1),
#              arrangeGrob(arrangeGrob(s.mouse.main, s.mouse.sub, nrow=1), arrangeGrob(s.chipmunk.main, s.chipmunk.sub, nrow=1), s.flyingsquirrel.main, nrow=1),
#              arrangeGrob(r2.mouse, r2.chipmunk, r2.flyingsquirrel, nrow=1),
#              arrangeGrob(rmse.mouse, rmse.chipmunk, rmse.flyingsquirrel, nrow=1),
#              nrow=4, align="v")
#              
# plot_grid(all.mouse, plot_grid(s.mouse.main, s.mouse.sub, nrow=1), r2.mouse, rmse.mouse, ncol=1, align="v", rel_heights = c(2,1,1,1.25))
# 
# grid.arrange(arrangeGrob(all.mouse, all.chipmunk, all.flyingsquirrel, nrow=1),
#           arrangeGrob(arrangeGrob(s.mouse.main, s.mouse.sub, nrow=1), arrangeGrob(s.chipmunk.main, s.chipmunk.sub, nrow=1), s.flyingsquirrel.main, nrow=1),
#           arrangeGrob(r2.mouse, r2.chipmunk, r2.flyingsquirrel, nrow=1),
#           arrangeGrob(rmse.mouse, rmse.chipmunk, rmse.flyingsquirrel, nrow=1),
#           nrow=4, align="v")

abs <- plot_grid(all.mouse, all.chipmunk, all.flyingsquirrel,
          s.mouse, s.chipmunk, s.flyingsquirrel.main, #s.flyingsquirrel, 
          r2.mouse, r2.chipmunk, r2.flyingsquirrel,
          rmse.mouse, rmse.chipmunk, rmse.flyingsquirrel,
          nrow=4, ncol=3,
          align="hv", rel_heights = c(2.25,1.8,0.9,1))
ggsave(abs, filename="Figures/nmix_fig_sep.tiff", height=15, width=15, units="in", dpi=400, compression="lzw")


#####################
#plot number of models that converged for each species and model

ggplot(data=rbind(mod.coef.all, data.frame(Estimate=0, Std..Error=NA, t.value=NA, Pr...t..=NA, sp="chipmunk", t=0, mod="base", r2=NA, rmse=NA, n=0)),
       aes(x=t, y=n, fill=mod)) +geom_hline(aes(yintercept=8), col="grey50",lty="dashed") + geom_bar(stat='identity', position="dodge") + 
  scale_fill_manual(values=c("#440154FF","#21908CFF","#FDE725FF"), name="model") +
  ylab("number sites with\nmodel convergence") + xlab("consolidation time (min)") +
  theme_bw(base_size=16) + theme(panel.grid = element_blank()) + facet_wrap(~sp)
ggsave(filename="Figures/nmix_fig_nconverged.tiff", height=4, width=8, dpi=400, compression="lzw")
