options(tidyverse.quiet = TRUE)
library(tidyverse)
library(DynSysSimR)

dataDirPath <- "../DynSysSimData"
evalDirPath <- file.path(dataDirPath, "evaluation")
forecastDirPath <- file.path(dataDirPath, "forecast")
truthDirPath <- file.path(dataDirPath, "truth")


# default parameters
parentSeed <- 1
truthName <- "L63_default"
noiseName <- "Gauss"
noiseScale <- 0
testDuration <- 20
nInSample <- 100
nOutOfSample <- 100
nLong <- 1e5
oosReps <- 100
nTrain <- 2^10
stepRate <- 2^5

truth <- readInfoAndData(truthDirPath, truthName)

nTest <- testDuration / truth$info$timeStep / stepRate
truthMean <- colMeans(truth$data[, -1])
truthSd <- sqrt(sum(colMeans((truth$data[, -1] - rep(truthMean, each=nrow(truth$data)))^2)))

cat(sprintf("randomSeed: %d...\n", parentSeed))
set.seed(parentSeed)
seeds <- sample.int(.Machine$integer.max, 8L)
names(seeds) <- c("truth", "oos", "noise", "fit", "assimilate", "predict", "longStart", "longCompare")

truthSample <- sampleTruth(truth$data, stepRate, nTrain, nTest, seed=seeds["truth"])

xOos <- sampleNTruthX(truth$data, stepRate, nOutOfSample+1, oosReps, seeds["oos"])

xTrainTruth <- truthSample$xTrain
xTest <- truthSample$xTest
tTest <- truthSample$tTest
tTrain <- truthSample$tTrain

cat(sprintf("noiseScale: %g...\n", noiseScale))
noiseSd <- noiseScale * truthSd
xTrainTrain <- observe(xTrainTruth, noiseSd, type=noiseName, seed=seeds["noise"])

nDeg <- 3

wsFull <- getWeightScheduleFull(nDeg = nDeg, d = d, kmax = 1)
pt <- proc.time()
model <- fitFullLoss(xTrainTrain, nDeg=nDeg, weightSchedule=wsFull, normalizationType="full", normalizationScale=0.2)
cat(sprintf("took %.2fs.\n", (proc.time()-pt)[3]))


 # Assimilation Error
assimilation <- drop(model$analysis)
stopifnot(nrow(assimilation) == nTrain)
assimilationErr <- sqrt(rowSums((xTrainTruth - assimilation)^2))
cat(sprintf("Assimilation Error: %.3f\n", sqrt(mean(assimilationErr^2))))

# Forecast Error
initialCondAssi <- assimilation[nTrain, , drop=FALSE]
predictionAssi <- predictFullLoss(model, initialCondAssi, nTest)
forcastAssiErr <- sqrt(rowSums((xTest - drop(predictionAssi))^2))
cat(sprintf("Forecast assi valid steps: %d\n", min(which(forcastAssiErr > 0.5*truthSd)-1)))

initialCondTrue <- xTrainTruth[nrow(xTrainTruth), , drop=FALSE]
predictionTrue <- predictFullLoss(model, initialCondTrue, nTest)
forcastTrueErr <- sqrt(rowSums((xTest - drop(predictionTrue))^2))
cat(sprintf("Forecast true valid steps: %d\n", min(which(forcastTrueErr > 0.5*truthSd)-1)))


wsCoef <- getWeightScheduleCoef(nDeg = nDeg, d = d, kmax = 10, type = "point")
pt <- proc.time()
model <- fitCoefLoss(xTrainTrain, nDeg=nDeg, weightSchedule=wsCoef, normalizationType="full", normalizationScale=0.2, verbose = TRUE)
cat(sprintf("took %.2fs.\n", (proc.time()-pt)[3]))


 # Assimilation Error
assimilation <- drop(model$analysis)
stopifnot(nrow(assimilation) == nTrain)
assimilationErr <- sqrt(rowSums((xTrainTruth - assimilation)^2))
cat(sprintf("Assimilation Error: %.3f\n", sqrt(mean(assimilationErr^2))))

# Forecast Error
initialCondAssi <- assimilation[nTrain, , drop=FALSE]
predictionAssi <- predictFullLoss(model, initialCondAssi, nTest)
forcastAssiErr <- sqrt(rowSums((xTest - drop(predictionAssi))^2))
cat(sprintf("Forecast assi valid steps: %d\n", min(which(forcastAssiErr > 0.5*truthSd)-1)))

initialCondTrue <- xTrainTruth[nrow(xTrainTruth), , drop=FALSE]
predictionTrue <- predictFullLoss(model, initialCondTrue, nTest)
forcastTrueErr <- sqrt(rowSums((xTest - drop(predictionTrue))^2))
cat(sprintf("Forecast true valid steps: %d\n", min(which(forcastTrueErr > 0.5*truthSd)-1)))
