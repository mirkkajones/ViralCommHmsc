setwd() 

library(Hmsc)

#Load unfitted Hmsc model created by script 1.
load(file = "unfitted_models.Rdata")

nChains = 4
samples = 250
thin = 5000

n = 1

m = models[[n]]
m = sampleMcmc(m, samples = samples, thin=thin,
               adaptNf=rep(ceiling(0.4*samples*thin),m$nr),
               transient = ceiling(0.5*samples*thin),
               nChains = nChains, nParallel = nChains)
models[[n]] = m

filename_out = paste("models_thin_", as.character(thin),
                     "_samples_", as.character(samples),
                     "_chains_",as.character(nChains),
                     ".Rdata",sep = "")

save(models, modelnames, file=filename_out)
