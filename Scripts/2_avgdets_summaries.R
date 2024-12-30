#########################
#########################
#data represents number of detections per camera in sherman grid for PEMA and tomahawk grid in TATO and GLOR
#SCR random on the y axis
#avg detections on the x axis

#########################
#load packages
require(ggplot2)
require(ggpubr)

#########################
#load data

# random <- read.table(file="Data_linearregs/CR_estimates.txt", sep=",", header=T) #incorrect PEMA grid numbers; grid 6,7,8 are supposed to be grid 7,8,9 
random <- read.table(file="Data_linearregs/CR_estimates_correctgridnums.txt", sep=",", header=T) #correct PEMA grid numbers; grids 7,8,9 
random <- random[random$model == "SCR random",] #has absolute values and scaled values
random <- random[,c("sp","grid",
                    "dhat","mean.random","sd.random","X2.5..random","X97.5..random", #remove extra columns
                    "dhat.scaled","X2.5.r.scaled","X97.5.r.scaled")]
names(random) <- c("sp","grid", "dhat.random","mean.random","sd.random","X2.5..random","X97.5..random","dhat.random.scaled","X2.5.r.scaled","X97.5.r.scaled")
table(random$sp)

#simple index file
avgdets <- read.table("Output_cam/output_avgdets_2023-03-10.txt", sep=",")
avgdets <- avgdets[avgdets$grid != 0,] #remove practice grid if in data
table(avgdets$sp)

#conversion table
sp.convert <- data.frame(sp=c("chipmunk","flying\nsquirrel","flying\nsquirrel","mouse"), acro=c("TATO","GLSA","GLOR","PEMA"))
sp.convert$sp <- factor(sp.convert$sp, levels=c("mouse","chipmunk","flying\nsquirrel"))

#########################
#probably not necessary since it's an index, not absolute abundance estimate
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
#1a. calculate slopes, RMSE, and R2
avgdets <- merge(sp.convert, avgdets, all.y=T, by.x="acro", by.y="sp") #convert sp from GLOR, PEMA, TATO to mouse, flying squirrel, chipmunk
avgdets <- merge(avgdets[,c("sp","grid","mean","sd","ncams","t","d")], random, all.x=T, by=c("sp","grid")) #merge avgdets with SCR random dataset
avgdets$d <- as.factor(avgdets$d)

#need to plot points and slopes separately for each species, consolidation window, and number of days of data to include
table(avgdets[,c("d","t")])
#      t
# d    0 15 60 1440
# 1   24 24 24   24
# 2   24 24 24   24
# 3   24 24 24   24
# 4   24 24 24   24
# 5   24 24 24   24
# 6   24 24 24   24
# 7   24 24 24   24
# 8   24 24 24   24
# 9   24 24 24   24

avgdets.scaled <- NULL
mod.coef.all <- NULL
mod.coef.scaled <- NULL

for(t in c(0, 15, 60, 1440)) #t in minutes
{
  print(t)
  for(sp in unique(avgdets$sp))
  {
    sp.nobreaks <- gsub(sp, pattern="\n", replacement="")
    print(sp.nobreaks)
    for(d in unique(avgdets$d))
    {
      data <- avgdets[avgdets$sp == sp & avgdets$d == d & avgdets$t == t,]
      
      #scale variables
      data$mean.scaled <- scale(data$mean) #center and scale independent var
      data$X2.5.scaled <- (data$mean - data$sd - attr(data$mean.scaled, "scaled:center"))/attr(data$mean.scaled, "scaled:scale")
      data$X97.5.scaled <- (data$mean + data$sd - attr(data$mean.scaled, "scaled:center"))/attr(data$mean.scaled, "scaled:scale")
      
      avgdets.scaled <- rbind(avgdets.scaled, data)
      
      #####
      #linear regression with non-scaled values
      l <- lm(dhat.random ~ mean, data=data)
      
      #look at residuals using plot
      tiff(filename=paste("Figures/residuals_plots/avgdet_",sp.nobreaks,"_t", t, "_d", d,".tiff", sep=""), height=8, width=8, units="in", res=300, compression="lzw")
      par(mfrow=c(2,2))
      plot(l)
      # plot(l.sqrt)
      # plot(l.log)
      dev.off()
      
      mod.coef <- data.frame(coef(summary(l)))[2,]
      mod.coef$sp <- sp
      mod.coef$d <- d
      mod.coef$t <- t
      mod.coef$r2 <- summary(l)$adj.r.squared
      mod.coef$rmse <- sqrt(mean(l$residuals^2))
      mod.coef.all <- rbind(mod.coef.all, mod.coef)
      
      ####
      #scaled linear regression
      l <- lm(dhat.random.scaled ~ mean.scaled, data=data) #linear regression
      
      #look at residuals using plot (this is the same as the non-scaled one)
      # tiff(filename=paste("Figures/residuals_plots/avgdet_",sp.nobreaks,"_t", t, "_d", d,"_scaled.tiff", sep=""), height=8, width=8, units="in", res=300, compression="lzw")
      # par(mfrow=c(2,2))
      # plot(l)
      # dev.off()
      
      mod.coef.s <- data.frame(coef(summary(l)))[2,]
      mod.coef.s$sp <- sp
      mod.coef.s$d <- d
      mod.coef.s$t <- t
      mod.coef.s$r2 <- summary(l)$adj.r.squared
      mod.coef.s$rmse <- sqrt(mean(l$residuals^2))
      mod.coef.scaled <- rbind(mod.coef.scaled, mod.coef.s)
      
    }
  }
}

# write.table(avgdets.scaled, file="Data_linearregs/avgdets_estimates.txt", sep=",", row.names=F)

######
#1b. plot un-scaled values together
# t <- 15
# sp <- "mouse"

mod.coef.all$sp <- factor(mod.coef.all$sp, levels=c("mouse","chipmunk","flyingsquirrel","flying\nsquirrel"))
mod.coef.all$t <- factor(mod.coef.all$t, levels=c(0,15,60,1440))

# for(t in c(0, 15, 60, 1440)) #t in minutes
# {
#   print(t)
  # avgdet.t <- avgdets.scaled[avgdets.scaled$sp == sp & avgdets.scaled$t == t,]
  # avgdet.t <- avgdets.scaled[avgdets.scaled$t == t,]
  # mod.coef.t <- mod.coef.all[mod.coef.all$t==t,]
  
  all.0 <- ggplot(data=avgdets.scaled[avgdets.scaled$t == 0,], aes(x=mean, y=dhat.random, col=d, group=d)) +
    geom_smooth(method="lm", aes(fill=d, lty=d), alpha=0.1) + #add ribbons
    geom_errorbar(aes(ymin=X2.5..random, ymax=X97.5..random), width=0.01, lty="dashed", linewidth=0.1) + #vertical error bars
    geom_errorbar(aes(xmin=mean-sd, xmax=mean+sd), width=0.01, lty="dashed", linewidth=0.1) +
    geom_point(size=3, pch=1, stroke=2) + #, alpha=0.5) +
    # scale_color_manual(values=c.vals) + scale_fill_manual(values=c.vals) +
    scale_color_brewer(palette="BrBG", name="days active") +
    scale_fill_brewer(palette="BrBG", name="days active") +
    scale_linetype(name="days active") +
    xlab("model estimate (detections/camera)") + ylab("SCR random (animals/ha)") +
    theme_bw(base_size=20) + theme(panel.grid = element_blank(), legend.position="top", strip.background = element_blank(), strip.text.x = element_blank()) + 
    facet_wrap(~sp, scales="free")
  all.0
  
  #slope
  s <- ggplot(mod.coef.all, aes(x=t, y=Estimate, col=d, shape=t)) + 
    geom_point(size=5, position=position_dodge(width=1)) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=1)) +
    # scale_color_manual(values=c.vals) +
    scale_color_brewer(palette="BrBG", name="days active") +
    geom_hline(yintercept=0, lty="dashed") +
    geom_hline(yintercept=1, lty="dashed") + 
    xlab("") + ylab("slope") + #scale_y_continuous(limits=c(0,NA)) +
    theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + 
    facet_wrap(~sp, scale="free_y")
  
  #r2
  r2 <- ggplot(mod.coef.all, aes(x=t, y=r2, col=d, shape=t)) + geom_point(size=5, position=position_dodge(width=1)) + 
    # scale_color_manual(values=c.vals) +
    scale_color_brewer(palette="BrBG", name="days active") +
    geom_hline(yintercept = 0, lty="dashed") + geom_hline(yintercept=1, lty="dashed") + 
    xlab("") + ylab("R2") +
    theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + 
    facet_wrap(~sp, scale="free_y")
  
  #rmse
  rmse <- ggplot(mod.coef.all, aes(x=t, y=rmse, col=d, shape=t)) + geom_point(size=5, position=position_dodge(width=1)) + 
    # scale_color_manual(values=c.vals) +
    scale_color_brewer(palette="BrBG", name="days active") +
    geom_hline(yintercept = 0, lty="dashed") + #geom_hline(yintercept=1, lty="dashed") + 
    xlab("model") + ylab("RMSE") +
    theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5)) + facet_wrap(~sp, scale="free_y")
  
  p <- ggarrange(all.0, s, r2, rmse, ncol=1, align="v", heights=c(0.5,0.25,0.15,0.25))
  p
  # ggsave(p, filename="Figures/avgdets_fig.tiff", height=15, width=15, units="in", dpi=400, compression="lzw")
  
# }

#########################
#1c. scale values and plot and calculate slopes, RMSE, and R2
  mod.coef.scaled$sp <- factor(mod.coef.scaled$sp, levels=c("mouse","chipmunk","flyingsquirrel","flying\nsquirrel"))
  mod.coef.scaled$t <- factor(mod.coef.scaled$t, levels=c(0,15,60,1440))
  
  all.0 <- ggplot(data=avgdets.scaled[avgdets.scaled$t == 0,], aes(x=mean.scaled, y=dhat.random, col=d, group=d)) +
    geom_smooth(method="lm", aes(fill=d, lty=d), alpha=0.1) + #add ribbons
    geom_errorbar(aes(ymin=X2.5..random, ymax=X97.5..random), width=0.01, lty="dashed", linewidth=0.1) + #vertical error bars
    geom_errorbar(aes(xmin=X2.5.scaled, xmax=X97.5.scaled), width=0.01, lty="dashed", linewidth=0.1) +
    geom_point(size=3, pch=1, stroke=2) + #, alpha=0.5) +
    # scale_color_manual(values=c.vals) + scale_fill_manual(values=c.vals) +
    scale_color_brewer(palette="BrBG", name="days active") +
    scale_fill_brewer(palette="BrBG", name="days active") +
    scale_linetype(name="days active") +
    xlab("scaled model estimate (detections/camera)") + ylab("SCR random (animals/ha)") +
    theme_bw(base_size=20) + theme(panel.grid = element_blank(), legend.position="top", strip.background = element_blank(), strip.text.x = element_blank()) + 
    facet_wrap(~sp, scales="free")
  all.0
  
  #slope
  s <- ggplot(mod.coef.scaled, aes(x=t, y=Estimate, col=d, shape=t)) + 
    geom_point(size=5, position=position_dodge(width=0.75)) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=0.75)) +
    # scale_color_manual(values=c.vals) +
    scale_color_brewer(palette="BrBG", name="days active") +
    geom_hline(yintercept=0, lty="dashed") +
    geom_hline(yintercept=1, lty="dashed") + 
    xlab("") + ylab("slope") + #scale_y_continuous(limits=c(0,NA)) +
    theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + 
    facet_wrap(~sp, scale="free_y")
  
  #r2
  r2 <- ggplot(mod.coef.scaled, aes(x=t, y=r2, col=d, shape=t)) + geom_point(size=5, position=position_dodge(width=0.75)) + 
    # scale_color_manual(values=c.vals) +
    scale_color_brewer(palette="BrBG", name="days active") +
    geom_hline(yintercept = 0, lty="dashed") + geom_hline(yintercept=1, lty="dashed") + 
    xlab("") + ylab("R2") +
    theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank()) + 
    facet_wrap(~sp, scale="free_y")
  
  #rmse
  rmse <- ggplot(mod.coef.scaled, aes(x=t, y=rmse, col=d, shape=t)) + geom_point(size=5, position=position_dodge(width=0.75)) + 
    # scale_color_manual(values=c.vals) +
    scale_color_brewer(palette="BrBG", name="days active") +
    geom_hline(yintercept = 0, lty="dashed") + #geom_hline(yintercept=1, lty="dashed") + 
    xlab("model") + ylab("RMSE") +
    theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5)) + facet_wrap(~sp, scale="free_y")
  
  p <- ggarrange(all.0, s, r2, rmse, ncol=1, align="v", heights=c(0.5,0.25,0.15,0.25))
  p
  # ggsave(p, filename="Figures/avgdets_scaled_fig.tiff", height=15, width=15, units="in", dpi=400, compression="lzw")
  
############
#2
  
# write.csv(mod.coef.scaled, file="Data_linearregs/avgdets_slopes.csv")

mod.coef.scaled %>% group_by(sp) %>%
  summarize(min(Estimate), max(Estimate), max(r2), min(rmse)) %>% mutate_if(is.numeric, round, digits=2)

#   sp                 `min(Estimate)` `max(Estimate)` `max(r2)` `min(rmse)`
#   <chr>                        <dbl>           <dbl>     <dbl>       <dbl>
# 1 "chipmunk"                    0.57            0.86      0.87        0.29
# 2 "flying\nsquirrel"            0.07            0.51      0.19        0.72
# 3 "mouse"                      -0.01            0.81      0.58        0.57
  
mod.coef.scaled[order(mod.coef.scaled$sp, mod.coef.scaled$t,mod.coef.scaled$d),] %>% mutate_if(is.numeric, round, 3)
  
  #                Estimate Std..Error t.value Pr...t..               sp d    t     r2  rmse
  # mean.scaled1      0.756      0.276   2.740    0.041            mouse 1    0  0.520 0.606
  # mean.scaled6      0.775      0.273   2.834    0.036            mouse 2    0  0.540 0.594
  # mean.scaled3      0.809      0.264   3.068    0.028            mouse 3    0  0.584 0.565
  # mean.scaled       0.809      0.278   2.904    0.034            mouse 4    0  0.553 0.585
  # mean.scaled7      0.773      0.307   2.519    0.053            mouse 5    0  0.471 0.637
  # mean.scaled4      0.736      0.324   2.269    0.073            mouse 6    0  0.409 0.673
  # mean.scaled8      0.717      0.330   2.170    0.082            mouse 7    0  0.382 0.688
  # mean.scaled2      0.702      0.335   2.097    0.090            mouse 8    0  0.361 0.700
  # mean.scaled5      0.694      0.339   2.049    0.096            mouse 9    0  0.348 0.707
  # mean.scaled28     0.720      0.286   2.516    0.053            mouse 1   15  0.470 0.637
  # mean.scaled33     0.734      0.289   2.541    0.052            mouse 2   15  0.476 0.634
  # mean.scaled30     0.751      0.305   2.462    0.057            mouse 3   15  0.458 0.645
  # mean.scaled27     0.679      0.372   1.824    0.128            mouse 4   15  0.279 0.743
  # mean.scaled34     0.604      0.412   1.467    0.202            mouse 5   15  0.161 0.802
  # mean.scaled31     0.544      0.415   1.313    0.246            mouse 6   15  0.108 0.827
  # mean.scaled35     0.527      0.408   1.290    0.253            mouse 7   15  0.100 0.831
  # mean.scaled29     0.494      0.407   1.214    0.279            mouse 8   15  0.073 0.843
  # mean.scaled32     0.479      0.412   1.162    0.298            mouse 9   15  0.055 0.851
  # mean.scaled55    -0.011      0.480  -0.022    0.983            mouse 1   60 -0.200 0.959
  # mean.scaled60     0.378      0.395   0.956    0.383            mouse 2   60 -0.014 0.882
  # mean.scaled57     0.451      0.394   1.144    0.304            mouse 3   60  0.049 0.854
  # mean.scaled54     0.198      0.462   0.428    0.686            mouse 4   60 -0.157 0.942
  # mean.scaled61     0.167      0.471   0.355    0.737            mouse 5   60 -0.170 0.947
  # mean.scaled58     0.228      0.449   0.509    0.633            mouse 6   60 -0.141 0.935
  # mean.scaled62     0.234      0.435   0.538    0.614            mouse 7   60 -0.134 0.933
  # mean.scaled56     0.195      0.430   0.453    0.669            mouse 8   60 -0.153 0.940
  # mean.scaled59     0.159      0.435   0.366    0.729            mouse 9   60 -0.169 0.947
  # mean.scaled82     0.471      0.382   1.234    0.272            mouse 1 1440  0.080 0.840
  # mean.scaled87     0.413      0.405   1.022    0.354            mouse 2 1440  0.007 0.873
  # mean.scaled84     0.437      0.382   1.144    0.304            mouse 3 1440  0.049 0.854
  # mean.scaled81     0.300      0.413   0.727    0.500            mouse 4 1440 -0.085 0.912
  # mean.scaled88     0.256      0.425   0.604    0.572            mouse 5 1440 -0.118 0.926
  # mean.scaled85     0.209      0.422   0.494    0.642            mouse 6 1440 -0.144 0.937
  # mean.scaled89     0.195      0.420   0.465    0.661            mouse 7 1440 -0.150 0.939
  # mean.scaled83     0.128      0.425   0.302    0.775            mouse 8 1440 -0.179 0.951
  # mean.scaled86     0.065      0.428   0.151    0.886            mouse 9 1440 -0.195 0.957
  
  # mean.scaled10     0.858      0.126   6.803    0.000         chipmunk 1    0  0.866 0.289
  # mean.scaled15     0.751      0.211   3.551    0.012         chipmunk 2    0  0.624 0.485
  # mean.scaled12     0.675      0.250   2.697    0.036         chipmunk 3    0  0.473 0.574
  # mean.scaled9      0.658      0.258   2.555    0.043         chipmunk 4    0  0.441 0.591
  # mean.scaled16     0.662      0.256   2.582    0.042         chipmunk 5    0  0.447 0.587
  # mean.scaled13     0.673      0.251   2.681    0.036         chipmunk 6    0  0.469 0.576
  # mean.scaled17     0.687      0.245   2.805    0.031         chipmunk 7    0  0.495 0.561
  # mean.scaled11     0.694      0.241   2.875    0.028         chipmunk 8    0  0.509 0.553
  # mean.scaled14     0.694      0.241   2.876    0.028         chipmunk 9    0  0.509 0.553
  # mean.scaled37     0.762      0.205   3.727    0.010         chipmunk 1   15  0.648 0.469
  # mean.scaled42     0.735      0.221   3.331    0.016         chipmunk 2   15  0.591 0.505
  # mean.scaled39     0.624      0.272   2.299    0.061         chipmunk 3   15  0.380 0.622
  # mean.scaled36     0.598      0.281   2.127    0.077         chipmunk 4   15  0.335 0.644
  # mean.scaled43     0.624      0.272   2.295    0.062         chipmunk 5   15  0.379 0.623
  # mean.scaled40     0.662      0.256   2.582    0.042         chipmunk 6   15  0.447 0.587
  # mean.scaled44     0.701      0.238   2.945    0.026         chipmunk 7   15  0.523 0.546
  # mean.scaled38     0.720      0.229   3.152    0.020         chipmunk 8   15  0.561 0.524
  # mean.scaled41     0.717      0.230   3.115    0.021         chipmunk 9   15  0.554 0.527
  # mean.scaled64     0.566      0.292   1.935    0.101         chipmunk 1   60  0.282 0.669
  # mean.scaled69     0.648      0.262   2.476    0.048         chipmunk 2   60  0.423 0.600
  # mean.scaled66     0.621      0.273   2.273    0.063         chipmunk 3   60  0.373 0.625
  # mean.scaled63     0.630      0.269   2.342    0.058         chipmunk 4   60  0.391 0.617
  # mean.scaled70     0.691      0.243   2.843    0.029         chipmunk 5   60  0.503 0.557
  # mean.scaled67     0.738      0.219   3.373    0.015         chipmunk 6   60  0.597 0.501
  # mean.scaled71     0.779      0.194   4.026    0.007         chipmunk 7   60  0.685 0.443
  # mean.scaled65     0.796      0.182   4.369    0.005         chipmunk 8   60  0.721 0.417
  # mean.scaled68     0.782      0.192   4.073    0.007         chipmunk 9   60  0.690 0.440
  # mean.scaled91     0.682      0.247   2.758    0.033         chipmunk 1 1440  0.486 0.567
  # mean.scaled96     0.724      0.226   3.200    0.019         chipmunk 2 1440  0.569 0.519
  # mean.scaled93     0.651      0.261   2.496    0.047         chipmunk 3 1440  0.428 0.598
  # mean.scaled90     0.628      0.270   2.329    0.059         chipmunk 4 1440  0.387 0.618
  # mean.scaled97     0.696      0.241   2.891    0.028         chipmunk 5 1440  0.513 0.552
  # mean.scaled94     0.752      0.211   3.565    0.012         chipmunk 6 1440  0.626 0.483
  # mean.scaled98     0.783      0.191   4.100    0.006         chipmunk 7 1440  0.693 0.438
  # mean.scaled92     0.800      0.179   4.478    0.004         chipmunk 8 1440  0.731 0.409
  # mean.scaled95     0.775      0.196   3.949    0.008         chipmunk 9 1440  0.676 0.450
  
  # mean.scaled19     0.507      0.315   1.610    0.159 flying\nsquirrel 1    0  0.185 0.722
  # mean.scaled24     0.441      0.331   1.331    0.232 flying\nsquirrel 2    0  0.099 0.759
  # mean.scaled21     0.400      0.340   1.178    0.283 flying\nsquirrel 3    0  0.053 0.778
  # mean.scaled18     0.368      0.346   1.065    0.328 flying\nsquirrel 4    0  0.019 0.792
  # mean.scaled25     0.369      0.345   1.067    0.327 flying\nsquirrel 5    0  0.019 0.792
  # mean.scaled22     0.369      0.345   1.067    0.327 flying\nsquirrel 6    0  0.019 0.792
  # mean.scaled26     0.376      0.344   1.091    0.317 flying\nsquirrel 7    0  0.027 0.789
  # mean.scaled20     0.389      0.342   1.140    0.298 flying\nsquirrel 8    0  0.041 0.783
  # mean.scaled23     0.384      0.343   1.122    0.305 flying\nsquirrel 9    0  0.036 0.785
  # mean.scaled46     0.364      0.346   1.051    0.334 flying\nsquirrel 1   15  0.015 0.793
  # mean.scaled51     0.277      0.359   0.770    0.471 flying\nsquirrel 2   15 -0.062 0.824
  # mean.scaled48     0.197      0.368   0.536    0.611 flying\nsquirrel 3   15 -0.113 0.843
  # mean.scaled45     0.136      0.373   0.364    0.728 flying\nsquirrel 4   15 -0.141 0.854
  # mean.scaled52     0.153      0.372   0.411    0.696 flying\nsquirrel 5   15 -0.135 0.852
  # mean.scaled49     0.157      0.371   0.422    0.687 flying\nsquirrel 6   15 -0.133 0.851
  # mean.scaled53     0.189      0.369   0.513    0.626 flying\nsquirrel 7   15 -0.118 0.845
  # mean.scaled47     0.242      0.364   0.664    0.531 flying\nsquirrel 8   15 -0.087 0.833
  # mean.scaled50     0.223      0.366   0.611    0.564 flying\nsquirrel 9   15 -0.098 0.838
  # mean.scaled73     0.328      0.352   0.932    0.387 flying\nsquirrel 1   60 -0.019 0.807
  # mean.scaled78     0.210      0.367   0.572    0.588 flying\nsquirrel 2   60 -0.106 0.841
  # mean.scaled75     0.135      0.373   0.362    0.730 flying\nsquirrel 3   60 -0.142 0.854
  # mean.scaled72     0.068      0.376   0.182    0.862 flying\nsquirrel 4   60 -0.160 0.861
  # mean.scaled79     0.094      0.375   0.251    0.810 flying\nsquirrel 5   60 -0.155 0.859
  # mean.scaled76     0.108      0.374   0.289    0.782 flying\nsquirrel 6   60 -0.151 0.857
  # mean.scaled80     0.149      0.372   0.400    0.703 flying\nsquirrel 7   60 -0.136 0.852
  # mean.scaled74     0.208      0.367   0.567    0.591 flying\nsquirrel 8   60 -0.107 0.841
  # mean.scaled77     0.186      0.369   0.504    0.632 flying\nsquirrel 9   60 -0.119 0.846
  # mean.scaled100    0.375      0.344   1.087    0.319 flying\nsquirrel 1 1440  0.025 0.789
  # mean.scaled105    0.303      0.356   0.851    0.427 flying\nsquirrel 2 1440 -0.041 0.816
  # mean.scaled102    0.205      0.367   0.557    0.598 flying\nsquirrel 3 1440 -0.109 0.842
  # mean.scaled99     0.119      0.374   0.317    0.762 flying\nsquirrel 4 1440 -0.147 0.856
  # mean.scaled106    0.150      0.372   0.404    0.700 flying\nsquirrel 5 1440 -0.136 0.852
  # mean.scaled103    0.164      0.371   0.441    0.675 flying\nsquirrel 6 1440 -0.130 0.850
  # mean.scaled107    0.221      0.366   0.604    0.568 flying\nsquirrel 7 1440 -0.100 0.838
  # mean.scaled101    0.295      0.357   0.826    0.440 flying\nsquirrel 8 1440 -0.047 0.818
  # mean.scaled104    0.262      0.361   0.725    0.495 flying\nsquirrel 9 1440 -0.073 0.828



##################

#plot each species separately and then combine at the end
main <- c(0, 15, 60, 1440)
sub <- c(0,15,60)
for(sp in unique(mod.coef.all$sp))
{
  print(sp)
  assign(paste("all.", gsub(sp, pattern="\n", replacement=""), sep=""), 
         ggplot(data=avgdets[avgdets$sp == sp & avgdets$t == 0,], aes(x=mean, y=dhat.random, col=d, group=d)) +
           geom_smooth(method="lm", aes(fill=d), alpha=0.1) + #add ribbons #, lty=mod
           geom_errorbar(aes(ymin=X2.5..random, ymax=X97.5..random), width=0.01, lty="dashed", linewidth=0.1) + #vertical error bars
           geom_errorbar(aes(xmin=mean-sd, xmax=mean+sd), width=0.01, lty="dashed", linewidth=0.1) +
           geom_point(size=3, pch=1, stroke=2) + #, alpha=0.5) +
           # scale_color_manual(values=c.vals) + scale_fill_manual(values=c.vals) +
           scale_color_brewer(palette="BrBG", name="days active") +
           scale_fill_brewer(palette="BrBG", name="days active") +
           scale_linetype(name="days active") +
           theme_bw(base_size=20) + theme(panel.grid = element_blank(), legend.position="top", strip.background = element_blank(), strip.text.x = element_blank()) + 
           theme(axis.title=element_blank()) + #remove all axis titles
           # facet_wrap(~sp, scales="free") +
           coord_cartesian(ylim=c(0, NA)))
  
  #slope
  assign(paste("s.", gsub(sp, pattern="\n", replacement=""), ".main", sep=""),
         ggplot(mod.coef.all[mod.coef.all$t %in% main & mod.coef.all$sp == sp,], aes(x=t, y=Estimate, col=d, shape=t)) + 
           geom_point(size=5, position=position_dodge(width=1)) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=1)) +
           scale_color_brewer(palette="BrBG", name="days active") +
           geom_hline(yintercept=0, lty="dashed") +
           geom_hline(yintercept=1, lty="dashed", col="grey50") + 
           xlab("") + ylab("slope") + #scale_y_continuous(limits=c(0,NA)) +
           theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank(), axis.title=element_blank()) #+
         # facet_wrap(~sp, scale="free_y")
  )
  assign(paste("s.", gsub(sp, pattern="\n", replacement=""), ".sub", sep=""),
         ggplot(mod.coef.all[mod.coef.all$t %in% sub & mod.coef.all$sp == sp,], aes(x=t, y=Estimate, col=d, shape=t)) + 
           geom_point(size=5, position=position_dodge(width=1)) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=1)) +
           scale_color_brewer(palette="BrBG", name="days active") +
           geom_hline(yintercept=0, lty="dashed") +
           geom_hline(yintercept=1, lty="dashed", col="grey50") + 
           xlab("") + ylab("") + #scale_y_continuous(limits=c(0,NA)) +
           theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), plot.margin = unit(c(0, 0, 0, 0),"inches")) #+ #, axis.text.x = element_blank(), axis.title=element_blank()
         # facet_wrap(~sp, scale="free_y")
  )
  
  assign(paste("r2.", gsub(sp, pattern="\n", replacement=""), sep=""),
         ggplot(mod.coef.all[mod.coef.all$sp == sp,], aes(x=t, y=r2, col=d, shape=t)) + geom_point(size=5, position=position_dodge(width=1)) + 
           scale_color_brewer(palette="BrBG", name="days active") +
           geom_hline(yintercept = 0, lty="dashed") + geom_hline(yintercept=1, lty="dashed") + 
           xlab("") + ylab("R2") +
           theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_blank(), axis.title=element_blank())) #+ 
  # facet_wrap(~sp, scale="free_y"))
  assign(paste("rmse.", gsub(sp, pattern="\n", replacement=""), sep=""),
         ggplot(mod.coef.all[mod.coef.all$sp == sp,], aes(x=t, y=rmse, col=d, shape=t)) + geom_point(size=5, position=position_dodge(width=1)) + 
           scale_color_brewer(palette="BrBG", name="days active") +
           geom_hline(yintercept = 0, lty="dashed") + #geom_hline(yintercept=1, lty="dashed") + 
           xlab("model") + ylab("RMSE") +
           theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), axis.text.x = element_text(angle = 90, vjust=0.5), axis.title=element_blank())) #+ facet_wrap(~sp, scale="free_y"))
}

all <- ggarrange(all.mouse, all.chipmunk, all.flyingsquirrel, nrow=1, common.legend = T)
all <- annotate_figure(all, left=text_grob("SCR random (animals/ha)", rot=90, vjust=1, size=20), bottom=text_grob("N-mixture density estimate (animals/ha)", size=20))

sub <- c(0,15)
sp <- "mouse"
assign(paste("s.", gsub(sp, pattern="\n", replacement=""), ".sub", sep=""),
       ggplot(mod.coef.all[mod.coef.all$t %in% sub & mod.coef.all$sp == sp,], aes(x=t, y=Estimate, col=d, shape=t)) + 
         geom_point(size=5, position=position_dodge(width=1)) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=1)) +
         scale_color_brewer(palette="BrBG", name="days active") +
         geom_hline(yintercept=0, lty="dashed") +
         geom_hline(yintercept=1, lty="dashed", col="grey50") + 
         xlab("") + ylab("") + #scale_y_continuous(limits=c(0,NA)) +
         theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), plot.margin = unit(c(0, 0, 0, 0),"inches")) #+ #, axis.text.x = element_blank(), axis.title=element_blank()
       # facet_wrap(~sp, scale="free_y")
)

sub <- c(60)
s.mouse.sub2 <-
  ggplot(mod.coef.all[mod.coef.all$t %in% sub & mod.coef.all$sp == sp,], aes(x=t, y=Estimate, col=d, shape=t)) + 
  geom_point(size=5, position=position_dodge(width=1)) + geom_errorbar(aes(ymin=Estimate - Std..Error, ymax=Estimate + Std..Error), width=0.1, linewidth=1, position=position_dodge(width=1)) +
  scale_color_brewer(palette="BrBG", name="days active") +
  geom_hline(yintercept=0, lty="dashed") +
  geom_hline(yintercept=1, lty="dashed", col="grey50") + 
  xlab("") + ylab("") + #scale_y_continuous(limits=c(0,NA)) +
  theme_bw(base_size=20) + theme(panel.grid=element_blank(), legend.position="none", strip.background = element_blank(), strip.text.x = element_blank(), plot.margin = unit(c(0, 0, 0, 0),"inches"))

s.mouse <- s.mouse.main + 
  annotation_custom(ggplotGrob(s.mouse.sub+rremove("ylab")+rremove("xlab")), xmin=0.5, xmax=2.25, ymin=12, ymax=125) + 
  annotation_custom(ggplotGrob(s.mouse.sub2+rremove("ylab")+rremove("xlab")), xmin=2.25, xmax=3.4, ymin=12, ymax=125)
s.chipmunk <- s.chipmunk.main + annotation_custom(ggplotGrob(s.chipmunk.sub+rremove("ylab")+rremove("xlab")), xmin=0.5, xmax=3.4, ymin=2.5, ymax=12)
s.flyingsquirrel <- s.flyingsquirrel.main + annotation_custom(ggplotGrob(s.flyingsquirrel.sub), xmin=0.5, xmax=3.5, ymin=0.5, ymax=2.25)

#combine in grid
require(cowplot)

abs <- plot_grid(all.mouse, all.chipmunk, all.flyingsquirrel,
                 s.mouse, s.chipmunk, s.flyingsquirrel.main, #s.flyingsquirrel, 
                 r2.mouse, r2.chipmunk, r2.flyingsquirrel,
                 rmse.mouse, rmse.chipmunk, rmse.flyingsquirrel,
                 nrow=4, ncol=3,
                 align="hv", rel_heights = c(2.25,1.8,0.9,1))
# ggsave(abs, filename="Figures/avgdets_scaled_fig_sep.tiff", height=15, width=15, units="in", dpi=400, compression="lzw")
