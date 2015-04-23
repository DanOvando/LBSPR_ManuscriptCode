WD <- "E:/GitRepos/LBSPR_ManuscriptCode"
setwd(WD)

# Load LBSPR Package
# install.packages() # install package from GitHub 
library(LBSPR)

SimPars <- LoadSimPars(PathtoSimFile=file.path(WD, "ExampleData"), SimParFileName="SimData", SimParExt=".csv", ind=1)

#############################################################
# Optimise kslope value for length-structured LHR GTG Model #
#############################################################

# kslope determines the slope of the linear function that control relative M for each GTG.
# M for mean GTG is equal to M (or M/K in LHR) and decrease and increases for smaller and larger
# GTG respectively with a linear model with slope = kslope 

OptimiseFitness <- function(logkslope, SimPars, Function, Mpar=NULL) {
  Kslope <- exp(logkslope) #+ 0.000000001
  SimPars$kslope <- Kslope
  return(Function(SimPars=SimPars, Mpar=0.1)$ObjFun)
}  

# Optimise kslope 
RunoptLHM <- optimise(OptimiseFitness, interval=log(c(0.000001, 0.2)), SimPars=SimPars, Function= SimMod_LHR)

# Test for patterns between optimal kslope and other parameters
# Random draws of all Simulation Parameters
TestKslopeFun <- function(Function) {
  SimPars$NGTG <- 60 # model is insensitive to NGTG as.integer(runif(1, 20, 100)) 
  SimPars$Linf <- runif(1, 30, 1000)
  SimPars$CVLinf <- runif(1, 0.08, 0.2)
  relL50 <- runif(1, 0.3, 0.9)
  SimPars$L50 <- relL50 * SimPars$Linf 
  SimPars$L95 <- runif(1, SimPars$L50, SimPars$Linf)
  MatDelta <-  SimPars$L95 - SimPars$L50 
  SimPars$L95 <- SimPars$L50 + MatDelta
  SimPars$Mpow <- runif(1, 0, 0.2)
  Mpar <- runif(1, 0.05, 0.5)
  SimPars$MK <- runif(1, 0.3, 4)
  
  if (SimPars$NGTG > SimPars$Linf) SimPars$NGTG <- as.integer(SimPars$Linf)
  SimPars$Linc <- 5
  if (SimPars$Linf/5 < 30) SimPars$Linc <- 1
  
  SimPars$FM <- 0 
  kslope <- exp(optimise(OptimiseFitness, interval=log(c(0.000001, 0.2)), SimPars=SimPars, Function= Function, Mpar=Mpar)$minimum)
  
  SimPars$kslope <- kslope
  RunFunc <- Function(SimPars=SimPars, Mpar=Mpar)
  OptSuccess <- TRUE
  if (RunFunc$Pen > 0) OptSuccess <- FALSE # Check if M or M/K is negative for any GTG
  
  Output <- c(kslope, SimPars$NGTG, SimPars$Linf, SimPars$CVLinf, SimPars$L50/SimPars$Linf, MatDelta, RunFunc$Mpar, SimPars$Mpow, SimPars$MK, RunFunc$kpar, OptSuccess)
  return(Output)
}

 # Life History Ratio Model Test 
# Uncomment the code below to re-run the test
# N <- 2000
# SaveRun <- matrix(NA, nrow=N, ncol=9)
# for (X in 1:N) {
  # print(X)
  # SaveRun[X,] <- TestKslopeFun(SimMod_LHR)
# }

# Fails <- which(SaveRun[,9] != 1)
# if (length(Fails) > 0) Run <- Run[-Fails,]
# Run <- SaveRun
# dput(Run, file="RCode/SaveKslopeRun.dat")

Run <- dget("RCode/SaveKslopeRun.dat")
colnames(Run) <- c("kslope", "NGTG", "Linf", "CVLinf", "RelL50", "MatDelta", "Mpow", "MK", "OptSuccess")
Run <- as.data.frame(Run)
par(mfrow=c(2,3), cex.lab=1.5, bty="l", mar=c(5,5,2,2), oma=c(1,1,1,1))
# plot(Run$NGTG, Run$kslope, xlab="NGTG", ylab="kslope") # no need for this one - fixed now 
plot(Run$Linf, Run$kslope, xlab="Linf", ylab="kslope")
plot(Run$CVLinf, Run$kslope, xlab="CVLinf", ylab="kslope")
plot(Run$RelL50, Run$kslope, xlab="L50/Linf", ylab="kslope")
plot(Run$MatDelta, Run$kslope, xlab="MatDelta", ylab="kslope")
plot(Run$Mpow, Run$kslope, xlab="Mpow", ylab="kslope")
plot(Run$MK, Run$kslope, xlab="MK", ylab="kslope")

Run$LogKSlope <- log(Run$kslope)

# Fit GAM 
require(mgcv)
# Mod1 <- gam(LogKSlope ~ s(NGTG, bs="cs") + s(Linf, bs="cs") + s(CVLinf, bs="cs") + s(RelL50, bs="cs") + s(MatDelta, bs="cs") + s(Mpow, bs="cs")+ s(MK, bs="cs"), data=Run)
Mod1 <- gam(LogKSlope ~ s(Linf, bs="cs") + s(CVLinf, bs="cs") + s(RelL50, bs="cs") + s(MatDelta, bs="cs") + s(Mpow, bs="cs")+ s(MK, bs="cs"), data=Run)
summary(Mod1)

Mod2 <- gam(LogKSlope ~ s(Linf, bs="cs") + s(CVLinf, bs="cs") + s(RelL50, bs="cs")  + s(MatDelta, bs="cs") + s(Mpow, bs="cs")+ s(MK, bs="cs"), data=Run)
summary(Mod2)
par(mfrow=c(2,2))
gam.check(Mod2)
par(mfrow=c(3,2))
plot.gam(Mod2, se=TRUE)
hist(resid(Mod2))

plot(Mod2, pages=1, residuals=T, pch=19, cex=0.25, scheme=1, col='#FF8000', shade=T,shade.col='black')

Run$LinfMK <- log(Run$Linf) - log(Run$MK)
plot(Run$LogKSlope, Run$LinfMK)

Mod3 <- gam(LogKSlope ~ s(LinfMK, bs="cs") + s(CVLinf, bs="cs") + s(RelL50, bs="cs") + s(MatDelta, bs="cs") + s(Mpow, bs="cs"), data=Run)
summary(Mod3)
plot(Run$LinfMK, Run$LogKSlope)
par(mfrow=c(2,2))
gam.check(Mod3)
par(mfrow=c(3,2))
plot.gam(Mod3, se=TRUE)
hist(resid(Mod3))

plot(Mod3, pages=1, residuals=T, pch=19, cex=0.25, scheme=1, col='#FF8000', shade=T,shade.col='black')

# Try Linear Model 
Mod4 <- gam(LogKSlope ~ CVLinf + RelL50 + Mpow + LinfMK, data=Run)
summary(Mod4)
anova(Mod3, Mod4, test="F")
AIC(Mod3, Mod4) # GAM is technically better but prefer a linear model if possible

Mod5 <- lm(LogKSlope ~ LinfMK + CVLinf + RelL50 + MatDelta + Mpow, data=Run)
summary(Mod5)
par(mfrow=c(2,3))
plot(Mod5)
hist(Mod5$residuals)

# Stick with Mod5 - Linear Model
KSlopeMod <- Mod5
save(KSlopeMod, file="RCode/KSlopeMod.rda") # Save Linear Model  
# Copy KSlopeMod over to the LBSPR Package data directory

# Test model prediction for a few examples
par(mfrow=c(5,5))
for (X in 1:25) {
  SimPars$NGTG <- 60 # as.integer(runif(1, 20, 100))
  SimPars$Linf <- runif(1, 30, 1000)
  SimPars$CVLinf <- runif(1, 0.08, 0.2)
  relL50 <- runif(1, 0.3, 0.9)
  SimPars$L50 <- relL50 * SimPars$Linf 
  SimPars$L95 <- runif(1, SimPars$L50, SimPars$Linf)
  MatDelta <-  SimPars$L95 - SimPars$L50 
  SimPars$Mpow <- runif(1, 0, 0.2)
  SimPars$MK <- runif(1, 0.3, 4)
  
  RunoptLHM <- optimise(OptimiseFitness, interval=log(c(0.000001, 0.2)), SimPars=SimPars, Function=SimMod_LHR)$minimum
  kslope <- exp(RunoptLHM)
  
  LMPredict <- PredictKSlope(SimPars)
  SimPars$kslope <- kslope
  LHM <- SimMod_LHR(SimPars)
  SimPars$kslope <- LMPredict
  LHMpredict <- SimMod_LHR(SimPars)
  SimPars$kslope <- 0
  LHMnoslope <- SimMod_LHR(SimPars)
  plot(LHM$FitPR, ylim=c(0, max(LHM$FitPR)*1.5))
  lines(LHMpredict$FitPR, ylim=c(0, max(LHMpredict$FitPR)), col="blue")
  lines(LHMnoslope$FitPR, ylim=c(0, max(LHMnoslope$FitPR)), col="red")
  print("********************")
  print(paste(X, "of 25"))
  print(cbind(kslope, LMPredict))
  print(t(SimPars[1:10]))
  print("********************")
}

# Result:
# Prediction model does relatively well for most parameter combinations. Certainly better than setting klope=0
# Tried other models for kslope, rather than linear model, e.g., cubic polynomial, but all optimial models for 
# approximating equal fitness across GTG resulted in a linear model 

##################################################################
# Test if kslope makes much difference to size structure and SPR #
##################################################################

# This test should be repeated for the paper and test if the value of kslope makes a difference to estimating 
# SPR. Here just to check how much difference between predicted and optimised values for kslope
CompSPR <- matrix(NA, nrow=25, ncol=3)
par(mfrow=c(5,5))
for (X in 1:25) {
  SimPars$NGTG <- as.integer(runif(1, 20, 100))
  SimPars$Linf <- runif(1, 30, 1000)
  SimPars$CVLinf <- runif(1, 0.08, 0.2)
  relL50 <- runif(1, 0.3, 0.9)
  SimPars$L50 <- relL50 * SimPars$Linf 
  SimPars$L95 <- SimPars$L50 + 0.1 * SimPars$Linf 
  SimPars$Mpow <- runif(1, 0, 0.2)
  SimPars$MK <- runif(1, 0.3, 4)
  SimPars$FM <- 0 # no fishing mortality when optimising kslope
  kslope <- exp(optimise(OptimiseFitness, interval=log(c(0.000001, 0.2)), SimPars=SimPars, Function=SimMod_LHR)$minimum)
  LMPredict <- PredictKSlope(SimPars, ModelLoc=ModelLoc)
  
  SimPars$FM <- runif(1, 0, 3)
  relSL50 <- runif(1, 0.3, 0.9)
  SimPars$SL50 <- relSL50 * SimPars$Linf 
  SimPars$SL95 <- SimPars$SL50 + 0.1 * SimPars$Linf 
  
  LHM <- SimMod_LHR(SimPars, kslope)
  LHMpredict <- SimMod_LHR(SimPars, LMPredict)
  LHMnoslope <- SimMod_LHR(SimPars, kslope=0)
  
  CompSPR[X,] <- cbind(LHM$SPR, LHMpredict$SPR, LHMnoslope$SPR)

  plot(LHM$LenMids, LHM$ExpLenCatchFished)
  lines(LHMpredict$LenMids, LHMpredict$ExpLenCatchFished, col="blue")
  lines(LHMnoslope$LenMids, LHMnoslope$ExpLenCatchFished, col="red")
}
colnames(CompSPR) <- c("Opt", "Predict", "Zero")
CompSPR <- as.data.frame(CompSPR)
CompSPR <- CompSPR[order(CompSPR[,1]),]
par(mfrow=c(1,2))
plot(CompSPR$Opt)
points(CompSPR$Predict, col="blue")
points(CompSPR$Zero, col="red")
plot(CompSPR$Opt - CompSPR$Predict, col="blue", pch=19, ylim=c(-0.02, 0.02))
points(CompSPR$Opt - CompSPR$Zero, col="red", pch=19)

# Result:
# predicted value of kslope returns size structure and estimated SPR very close to optimised value 
# kslope=0 returns a close value of SPR but affects size structure, sometimes quite dramatically
# Recommend predicting kslope from linear model and testing for sensitivity with LBSPR model


