#########################
#########################

#SCR random on the y axis
#TTE and STE estimates on x axis

#########################
#load packages
require(dplyr)

require(ggplot2)
require(ggpubr)
require(viridis)

#########################
#load data

#conversion table
# sp.convert <- data.frame(sp=c("chipmunk","flying\nsquirrel","flying\nsquirrel","mouse"), acro=c("TATO","GLSA","GLOR","PEMA"))
sp.convert <- data.frame(sp=c("chipmunk","flying\nsquirrel","mouse"), sp2=c("TOWNSENDS CHIPMUNK","FLYING SQUIRREL","DEER MOUSE"))
sp.convert$sp <- factor(sp.convert$sp, levels=c("mouse","chipmunk","flying\nsquirrel"))

#random SCR estimates
# random <- read.table(file="Data_linearregs/CR_estimates.txt", sep=",", header=T) #wrong grid numbers for PEMA
random <- read.table(file="Data_linearregs/CR_estimates_correctgridnums.txt", sep=",", header=T)
random <- random[random$model == "SCR random",] #has absolute values and scaled values
random <- random[,c("sp","grid",
                    "dhat","mean.random","sd.random","X2.5..random","X97.5..random", #remove extra columns
                    "dhat.scaled","X2.5.r.scaled","X97.5.r.scaled")]
names(random) <- c("sp","grid", "dhat.random","mean.random","sd.random","X2.5..random","X97.5..random","dhat.random.scaled","X2.5.r.scaled","X97.5.r.scaled")
table(random$sp)

tte <- read.table("Output_cam/output_STE_TTE_model_results_2020-08-26.txt", sep=",", header=T)

#t <- params[params$sp == sp,]$t # duration of each sampling occasion for STE
#m <- params[params$sp == sp,]$m # animal speed in m per hour for TTE
#sa <- params[params$sp == sp,]$sa # study area size

####################

tte <- merge(sp.convert, tte, by.y="sp", by.x="sp2", all.y=T) #factor species names
tte$sp <- factor(tte$sp, levels=c("mouse","chipmunk","flying\nsquirrel"))
tte$t <- factor(tte$t, levels=c(1,5,15,30,60,1440)) #factor delta 
tte$model <- factor(tte$model, levels=c("tte","ste")) #factor model

tte <- merge(tte, random, by.x=c("g","sp"), by.y=c("grid","sp"))

#convert to density
tte$D <- tte$N/(tte$sa/10000)
tte$X2.5. <- tte$LCI/(tte$sa/10000)
tte$X97.5. <- tte$UCI/(tte$sa/10000)

#scale tte mean values
#only scale values that converged
tte$D.scaled <- NA
tte$X2.5.scaled <- NA
tte$X97.5.scaled <- NA

mod.coef.all <- data.frame()
mod.coef.scaled <- data.frame()

for(t in unique(tte$t))
{
  print(t)
  for(sp in unique(tte$sp))
  {
    print(sp)
    for(m in unique(tte$model))
    {
      print(m)
      if(nrow(tte[tte$t == t & tte$sp == sp & tte$model == m,]) > 0)
      {
        ########
        #not scaled
        #linear regression with non-scaled values
        l <- lm(dhat.random ~ D, data=tte[tte$t == t & tte$sp == sp & tte$model == m,])
        
        mod.coef <- data.frame(coef(summary(l)))[2,]
        mod.coef$sp <- sp
        mod.coef$t <- t
        mod.coef$mod <- m
        mod.coef$r2 <- summary(l)$adj.r.squared
        mod.coef$rmse <- sqrt(mean(l$residuals^2))
        mod.coef.all <- rbind(mod.coef.all, mod.coef)
        
        ########
        #scaled
        s <- scale(tte[tte$t == t & tte$sp == sp  & tte$model == m,]$D)
        tte[tte$t == t & tte$sp == sp & tte$model == m,]$D.scaled <- s
        
        tte[tte$t == t & tte$sp == sp & tte$model == m,]$X2.5.scaled <- (tte[tte$t == t & tte$sp == sp & tte$model == m,]$X2.5. - attr(s, "scaled:center"))/attr(s, "scaled:scale")
        tte[tte$t == t & tte$sp == sp & tte$model == m,]$X97.5.scaled <- (tte[tte$t == t & tte$sp == sp & tte$model == m,]$X97.5. - attr(s, "scaled:center"))/attr(s, "scaled:scale")
        
        if(nrow(tte[tte$t == t & tte$sp == sp & tte$model == m,]) == 1)
        {
          tte[tte$t == t & tte$sp == sp & tte$model == m,]$D.scaled <- 0
        }
        
        l <- lm(dhat.random.scaled ~ D.scaled, data=tte[tte$t == t & tte$sp == sp & tte$model == m,]) #linear regression
        
        mod.coef.s <- data.frame(coef(summary(l)))[2,]
        mod.coef.s$sp <- sp
        mod.coef.s$t <- t
        mod.coef.s$mod <- m
        mod.coef.s$r2 <- summary(l)$adj.r.squared
        mod.coef.s$rmse <- sqrt(mean(l$residuals^2))
        mod.coef.scaled <- rbind(mod.coef.scaled, mod.coef.s)
        
      }
    }
  }
}

# write.table(tte, file="Data_linearregs/tte_estimates.txt", sep=",", row.names=F)

######
#1b. plot unscaled values together

mod.coef.all$model <- paste(mod.coef.all$mod, mod.coef.all$t, sep=" ")
mod.coef.all$model <- factor(mod.coef.all$model, levels=c("tte 1440","ste 1","ste 5","ste 15","ste 30","ste 60"))
mod.coef.all$sp <- factor(mod.coef.all$sp, levels=c("mouse","chipmunk","flyingsquirrel","flying\nsquirrel"))
mod.coef.all$t <- factor(mod.coef.all$t, levels=c(1,5,15,30,60,1440))

tte$mod <- paste(tte$model, tte$t, sep=" ")
tte$mod <- factor(tte$mod, levels=c("tte 1440","ste 1","ste 5","ste 15","ste 30","ste 60"))

all.0 <- ggplot(data=tte, aes(x=D, y=dhat.random, col=mod, shape=mod)) +
  geom_smooth(method="lm", aes(fill=mod, lty=mod), alpha=0.1) + #add ribbons
  geom_errorbar(aes(ymin=X2.5..random, ymax=X97.5..random), width=0.01, lty="dashed", linewidth=0.1) + #vertical error bars
  geom_errorbar(aes(xmin=X2.5., xmax=X97.5.), width=0.01, lty="dashed", linewidth=0.1) +
  geom_point(size=3, stroke=2) + #, alpha=0.5) +
  # scale_color_manual(values=c.vals) + scale_fill_manual(values=c.vals) +
  scale_color_jco(name="model") + scale_fill_jco(name="model") +
  scale_linetype(name="model") +
  scale_shape_manual(name="model", values=c(2, rep(1, 5))) +
  xlab("model estimate (animals/ha)") + ylab("SCR random (animals/ha)") +
  theme_bw(base_size=20) + theme(panel.grid = element_blank(), legend.position="top", strip.background = element_blank(), strip.text.x = element_blank()) + 
  facet_wrap(~sp, scales="free") +
  coord_cartesian(ylim=c(0, NA))
# coord_cartesian(ylim=c(0, 150))
all.0

#slope
s <- ggplot(mod.coef.all, aes(x=model, y=Estimate, col=model, shape=mod)) + 
  geom_point(size=5, position=position_dodge(width=1)) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=1)) +
  scale_color_jco(name="model") +
  geom_hline(yintercept=0, lty="dashed") +
  geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("slope") + #scale_y_continuous(limits=c(0,NA)) +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + 
  facet_wrap(~sp, scale="free_y")

#r2
r2 <- ggplot(mod.coef.all, aes(x=model, y=r2, col=model, shape=mod)) + geom_point(size=5, position=position_dodge(width=1)) + 
  # scale_color_manual(values=c.vals) +
  scale_color_jco(name="model") +
  geom_hline(yintercept = 0, lty="dashed") + geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("R2") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + 
  facet_wrap(~sp, scale="free_y")

#rmse
rmse <- ggplot(mod.coef.all, aes(x=model, y=rmse, col=model, shape=mod)) + geom_point(size=5, position=position_dodge(width=1)) + 
  # scale_color_manual(values=c.vals) +
  scale_color_jco(name="model") +
  geom_hline(yintercept = 0, lty="dashed") + #geom_hline(yintercept=1, lty="dashed") + 
  xlab("model") + ylab("RMSE") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5)) + facet_wrap(~sp, scale="free_y")

p <- ggarrange(all.0, s, r2, rmse, ncol=1, align="v", heights=c(0.5,0.25,0.15,0.25))
p
# ggsave(p, filename="Figures/tte_fig.tiff", height=15, width=15, units="in", dpi=400, compression="lzw")

################

#########################
#1c. scale values and plot and calculate slopes, RMSE, and R2
mod.coef.scaled$sp <- factor(mod.coef.scaled$sp, levels=c("mouse","chipmunk","flyingsquirrel","flying\nsquirrel"))
mod.coef.scaled$t <- factor(mod.coef.scaled$t, levels=c(1,5,15,30,60,1440))

mod.coef.scaled$model <- paste(mod.coef.scaled$mod, mod.coef.scaled$t, sep=" ")
mod.coef.scaled$model <- factor(mod.coef.scaled$model, levels=c("tte 1440","ste 1","ste 5","ste 15","ste 30","ste 60"))


all.0 <- ggplot(data=tte, aes(x=D.scaled, y=dhat.random, col=mod, shape=mod)) +
  geom_smooth(data=tte, method="lm", aes(fill=mod, lty=mod), alpha=0.1) + #add ribbons
  geom_errorbar(aes(ymin=X2.5..random, ymax=X97.5..random), width=0.01, lty="dashed", linewidth=0.1) + #vertical error bars
  geom_errorbar(aes(xmin=X2.5.scaled, xmax=X97.5.scaled), width=0.01, lty="dashed", linewidth=0.1) +
  geom_point(size=3, stroke=1) + #, alpha=0.5) +
  scale_color_jco(name="model") + scale_fill_jco(name="model") +
  scale_linetype(name="model") +
  scale_shape_manual(name="model", values=c(2, rep(1, 5))) +
  xlab("scaled model estimate (detections/camera)") + ylab("SCR random (animals/ha)") +
  theme_bw(base_size=20) + theme(panel.grid = element_blank(), legend.position="top", strip.background = element_blank(), strip.text.x = element_blank()) + 
  facet_wrap(~sp, scales="free") +
  coord_cartesian(ylim=c(0, NA))
# coord_cartesian(ylim=c(0, 150))
all.0

#slope
s <- ggplot(mod.coef.scaled, aes(x=model, y=Estimate, col=model, shape=mod)) + 
  geom_point(size=5, position=position_dodge(width=0.75)) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=0.75)) +
  scale_color_jco(name="model") +
  geom_hline(yintercept=0, lty="dashed") +
  geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("slope") + #scale_y_continuous(limits=c(0,NA)) +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + 
  facet_wrap(~sp, scale="free_y")

#r2
r2 <- ggplot(mod.coef.scaled, aes(x=model, y=r2, col=model, shape=mod)) + geom_point(size=5, position=position_dodge(width=0.75)) + 
  scale_color_jco(name="model") +
  geom_hline(yintercept = 0, lty="dashed") + geom_hline(yintercept=1, lty="dashed") + 
  xlab("") + ylab("R2") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + 
  facet_wrap(~sp, scale="free_y")

#rmse
rmse <- ggplot(mod.coef.scaled, aes(x=model, y=rmse, col=model, shape=mod)) + geom_point(size=5, position=position_dodge(width=0.75)) + 
  scale_color_jco(name="model") +
  geom_hline(yintercept = 0, lty="dashed") + #geom_hline(yintercept=1, lty="dashed") + 
  xlab("model") + ylab("RMSE") +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5)) + facet_wrap(~sp, scale="free_y")

p <- ggarrange(all.0, s, r2, rmse, ncol=1, align="v", heights=c(0.5,0.25,0.15,0.25))
p
# ggsave(p, filename="Figures/tte_scaled_fig.tiff", height=15, width=15, units="in", dpi=400, compression="lzw")

write.csv(mod.coef.scaled, file="Data_linearregs/tte_slopes.csv")

mod.coef.scaled %>% group_by(sp) %>%
  summarize(min(Estimate, na.rm=T), max(Estimate, na.rm=T), max(r2), min(rmse))

#   sp                 `min(Estimate, na.rm = T)` `max(Estimate, na.rm = T)` `max(r2)` `min(rmse)`
#   <fct>                                   <dbl>                      <dbl>     <dbl>       <dbl>
# 1 "mouse"                                -0.684                      0.462    0.323        0.721
# 2 "chipmunk"                              0.454                      0.588    0.319        0.652
# 3 "flying\nsquirrel"                     -0.383                      0.438    0.0962       0.760

mod.coef.scaled[order(mod.coef.scaled$sp, mod.coef.scaled$mod, mod.coef.scaled$t),] %>% mutate_if(is.numeric, round, 3)

#            Estimate Std..Error t.value Pr...t..               sp    t mod     r2  rmse    model
# D.scaled8    -0.582      0.383  -1.517    0.190            mouse   15 ste  0.178 0.794   ste 15
# D.scaled14    0.262      0.448   0.584    0.585            mouse   60 ste -0.123 0.928   ste 60
# D.scaled2    -0.684      0.348  -1.964    0.107            mouse    1 ste  0.323 0.721    ste 1
# D.scaled5    -0.054      0.463  -0.116    0.912            mouse    5 ste -0.197 0.958    ste 5
# D.scaled11   -0.218      0.453  -0.481    0.651            mouse   30 ste -0.147 0.938   ste 30
# D.scaled17    0.462      0.415   1.115    0.316            mouse 1440 tte  0.039 0.859 tte 1440

# D.scaled6     0.518      0.306   1.690    0.142         chipmunk   15 ste  0.210 0.702   ste 15
# D.scaled12    0.578      0.288   2.009    0.091         chipmunk   60 ste  0.302 0.660   ste 60
# D.scaled      0.458      0.322   1.423    0.204         chipmunk    1 ste  0.128 0.738    ste 1
# D.scaled3     0.588      0.285   2.068    0.084         chipmunk    5 ste  0.319 0.652    ste 5
# D.scaled9     0.538      0.301   1.791    0.124         chipmunk   30 ste  0.240 0.689   ste 30
# D.scaled15    0.454      0.323   1.404    0.210         chipmunk 1440 tte  0.122 0.740 tte 1440

# D.scaled7     0.072      0.376   0.191    0.855 flying\nsquirrel   15 ste -0.160 0.861   ste 15
# D.scaled13    0.438      0.332   1.321    0.235 flying\nsquirrel   60 ste  0.096 0.760   ste 60
# D.scaled1    -0.065      0.491  -0.131    0.902 flying\nsquirrel    1 ste -0.245 0.897    ste 1
# D.scaled4    -0.074      0.376  -0.198    0.850 flying\nsquirrel    5 ste -0.159 0.861    ste 5
# D.scaled10    0.268      0.361   0.742    0.486 flying\nsquirrel   30 ste -0.069 0.826   ste 30
# D.scaled16   -0.383      0.343  -1.116    0.307 flying\nsquirrel 1440 tte  0.034 0.786 tte 1440
