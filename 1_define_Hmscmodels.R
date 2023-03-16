setwd() 

library(Hmsc)

load(file="allData.R")

# Exclude 6 plant individuals (data rows) that have missing values in those XData columns that define
# viral taxon presence/absence at early summer time point 1.
excludeplants1 = which(is.na(rowSums(X[, c("CLOSTERO_T1", "ENAMO_T1", "BETAPARTITI_T1", "PILV_T1", "CAULIMO_T1")])))

X = droplevels(X[-excludeplants1,])
Y = Y[-excludeplants1,]
S = data.frame(S[-excludeplants1,])

# Exclude 10 plant individuals with missing plant area data.
excludeplants2 = which(is.na(X$Plant.area))

X = droplevels(X[-excludeplants2,])
Y = Y[-excludeplants2,]
S = data.frame(S[-excludeplants2,])

# Define several columns in X as categorical (i.e., factors).
X$Genotype = as.factor(X$Genotype)
X$CLOSTERO_T1 = as.factor(X$CLOSTERO_T1)
X$BETAPARTITI_T1 = as.factor(X$BETAPARTITI_T1)
X$PILV_T1 = as.factor(X$PILV_T1)

names(X)[7:8] = c("HERBIVORY", "Plant.area")

# Define the Hmsc model Xformula: 
# X covariates are plant genotype, presence/absence of three viral taxa at time point 1, 
# herbivory and area of the host plant, and the number (if any) of samples used to define late summer 
# viral presence/absence that were missing (due e.g. to unsuccessful DNA extraction):
# Note: ENAMO and CAULIMO viruses are always absent at time 1 - hence are not included as X-covariates.

XFormula = ~ Genotype + CLOSTERO_T1 + BETAPARTITI_T1 + PILV_T1 + HERBIVORY + Plant.area + viral_NAs

# Define five XData matrices (one per modelled virus) as a list, because the number of missing sample values 
# varies among taxa.
XData = list()
for(n in 1:ncol(Y))
{
  Xdataprep = X
  if(n > 3)
  {Xdataprep$viral_NAs == 0}
  XData[[n]] = Xdataprep
}

# Define the StudyDesign matrix and two random effects: plant individual ID and plant population of origin:
studyDesign = data.frame(plant = as.factor(S$PlantID), population = as.factor(S$Population))

Population = studyDesign$population
rL.population = HmscRandomLevel(units = levels(Population))

Plant = studyDesign$plant
rL.plant = HmscRandomLevel(units = levels(Plant))

#Define Hmsc model structure
m1 = Hmsc(Y=Y, XData = XData,  XFormula = XFormula,
          distr="probit",
          studyDesign=studyDesign,
          ranLevels={list("plant" = rL.plant, "population" = rL.population)})

models = list(m1)
modelnames = c("presence_absence")

#Save the unfitted model structure and model name
save(models,modelnames,file = "unfitted_models.Rdata")
