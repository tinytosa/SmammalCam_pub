
################
#get minimum number known alive from SCR_allgrids data object

mnka.all <- NULL

sp.cr <- "PEMA"
sp.cr <- "GLSA"
sp.cr <- "TATO"

load(paste("RData/CR-SCR/SCR_allgrids_jagsargs_", sp.cr, ".RData", sep=""))

mnka <- NULL
for(g in 1:length(data$sf$caphist))
{
  mnka <- rbind(mnka, dim(data$sf$caphist[[g]]))
}

mnka <- data.frame(mnka)
names(mnka) <- c("mnka","ntraps","noccasions")
mnka$sp <- sp.cr
mnka$g <- 1:nrow(mnka)

mnka.all <- rbind(mnka, mnka.all)

#################
#do this after running through all the sp.cr
mnka.all$model <- "mnka"

write.table(mnka.all, file="Output/MNKA.txt", sep=",")
