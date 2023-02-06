setwd()

library(Hmsc)
library(vioplot)

load("models_thin_5000_samples_250_chains_4.Rdata")
m = models[[1]]

#Check beta parameter convergence
mpost = convertToCodaObject(m, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
summary(psrf.beta[,1])
vioplot(psrf.beta[,1],col="cyan",ylim=c(0.9,1.1),main="psrf(beta)")

#Check model explanatory performance (Tjur R2, AUC)
MF = evaluateModelFit(hM=m, predY=preds)
TjurR2 = MF$TjurR2
AUC = MF$AUC

# Run 10-fold cross-validation to assess model predictive performance (cv-Tjur R2, cv-AUC)
# Warning: Slow!
nChains = 4
samples = 250
thin = 5000

MF = list()
MFCV = list()
WAIC = list()

for(n in 1:1){
  m = models[[1]]
  preds = computePredictedValues(m)
  MF[[n]] = evaluateModelFit(hM=m, predY=preds)
  partition = createPartition(m, nfolds = 10)
  preds = computePredictedValues(m, partition=partition, nParallel = nChains)
  MFCV[[n]] = evaluateModelFit(hM=m, predY=preds)
  WAIC[[n]] = computeWAIC(m)       
}

filename_out = paste("MF_models_thin_", as.character(thin),
                     "_samples_", as.character(samples),
                     "_chains_",as.character(nChains),
                     ".Rdata",sep = "")

save(MF, MFCV, WAIC, modelnames, file = filename_out)

# Generate model fit summary Table 1:
cvTjurR2 = MFCV[[1]]$TjurR2
cvAUC = MFCV[[1]]$AUC
Table1 = data.frame(round(cbind(TjurR2, AUC, cvTjurR2, cvAUC),2))
rownames(Table1) = colnames(m$Y)
write.csv2(Table1, "Table1.csv")

