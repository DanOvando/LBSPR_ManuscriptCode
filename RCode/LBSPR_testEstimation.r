WD <- "E:/GitRepos/LBSPR_ManuscriptCode"
setwd(WD)

# Load LBSPR Package
# install.packages() # install package from GitHub 
library(LBSPR)

# This code tests the functionality of the ADMB estimation model, and checks that the simulation model and 
# estimation model return the same values when given the same inputs.

SimPars <- LoadSimPars(PathtoSimFile=file.path(WD, "ExampleData"), SimParFileName="SimData", SimParExt=".csv", ind=1)

AssessPars <- LoadAssessPars(PathtoAssessFile=file.path(WD, "ExampleData"), AssessParFileName="SimData", AssessParExt=".csv", ind=1)

####################################
# Test Model for set of parameters #
####################################
# SimPars and AssessPars are identical 
SimPars$FM <- 1
SimPars$MK <- AssessPars$MK <- 1.5
AssessPars$kslope <- SimPars$kslope <- PredictKSlope(SimPars)
# Generate Length Data 
LHRMod <- SimMod_LHR(SimPars)

# Run Assessmment Code 
LenFreq <- LHRMod$ExpLenCatchFished 
LenMids <- LHRMod$LenMids 
ADMBDir <- "E:/LBSPR_ADMB_Src/Copy1"
runMod <- RunLBSPRAssess(AssessPars, LenFreq, LenMids, ADMBDir, ExName="lbspr", showOutput=TRUE, MaxCount=5, ADMBRead=NULL)

plot(LHRMod$LenMids, LHRMod$ExpLenCatchFished)
lines(runMod$Bins, runMod$Pred)

cbind(Sims=c(SimPars$FM,LHRMod$SPR, SimPars$SL50, SimPars$SL95), Ests=runMod$Estimates[,2])

# Result:
# ADMB estimation model and R simulation model return exactly the same value... phew!

#################################
# Look at sensitivity to kslope #
#################################
# Initial test - more to follow
SimPars$MK <- AssessPars$MK <- 3
Multi <- 1.5
plot(c(0,Multi), c(0,1), type="n", xlab="Est kslope / Sim kslope", ylab="SPR", bty="l", main=paste("M/K = ", SimPars$MK))

FMVec <- seq(from=0, to=4, length=8) 
for (X in seq_along(FMVec)) {
  SimPars$FM <- FMVec[X]

  # Generate Length Data 
  LHRMod <- SimMod_LHR(SimPars)
  SimPars$kslope <- PredictKSlope(SimPars)
  
  LenFreq <- LHRMod$ExpLenCatchFished 
  LenMids <- LHRMod$LenMids 
  ADMBDir <- "E:/LBSPR_ADMB_Src/Copy1"
  
  N <- 10 
  KSlopeVec <- seq(from=0, to=SimPars$kslope*Multi, length=N)
  Saves <- matrix(NA, nrow=N, ncol=2)
  for (It in 1:N) {
    # Run Assessmment Code 
    AssessPars$kslope <- KSlopeVec[It]
    runMod <- RunLBSPRAssess(AssessPars, LenFreq, LenMids, ADMBDir, ExName="lbspr", showOutput=TRUE, MaxCount=5, ADMBRead=NULL)
    Saves[It,] <- runMod$Estimates[1:2,2]
  }
  
  lines(KSlopeVec/SimPars$kslope, Saves[,2], col=X)
}
abline(h=LHRMod$SPR, lty=2, lwd=2)
abline(v=1, lty=2, lwd=2)
