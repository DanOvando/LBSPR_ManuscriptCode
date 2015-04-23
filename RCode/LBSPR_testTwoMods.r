WD <- "E:/GitRepos/LBSPR_ManuscriptCode"
setwd(WD)

# Load LBSPR Package
install.packages() # install package from GitHub 
library(LBSPR)

SimPars <- LoadSimPars(PathtoSimFile=file.path(WD, "ExampleData"), SimParFileName="SimData", SimParExt=".csv", ind=1)

######################################################
# Test that the two models provide identical results #
######################################################
LHRMod <- SimMod_LHR(SimPars)
AgeMod <- SimMod_AgeStruc(SimPars, Mpar=0.1)

LHRMod$SPR
AgeMod$SPR # Should be identical or very close

LHRMod$Yield
AgeMod$Yield    
# Won't neccessarily be the same but should be same relative
# to maximum yield for each model 

# Test for maximum yield 
LHRYield <- sapply(seq(0, 2, length=100), function (X) {
  SimPars$FM <- X
  SimMod_LHR(SimPars)$Yield
})
AgeYield <- sapply(seq(0, 2, length=100), function (X) {
  SimPars$FM <- X
  SimMod_AgeStruc(SimPars, Mpar=0.1)$Yield
})

plot(seq(0, 2, length=100), LHRYield/max(LHRYield), xlab="F/M", ylab="Rel. Yield")
lines(seq(0, 2, length=100), AgeYield/max(AgeYield))

# Compare the two models - same size structure?
SimPars$FM <- 1 
LHRMod <- SimMod_LHR(SimPars)
AgeMod <- SimMod_AgeStruc(SimPars, Mpar=0.1)
plot(LHRMod$LenMids, LHRMod$ExpLenCatchFished)
lines(AgeMod$LenMids, AgeMod$ExpLenCatchFished)

