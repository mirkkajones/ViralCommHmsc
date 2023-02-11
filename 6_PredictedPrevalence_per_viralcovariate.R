setwd()

library(Hmsc)
library(ggplot2)
library(gridExtra)
library(abind)

load("models/models_thin_5000_samples_250_chains_4.Rdata")
m = models[[1]]

#set the random number seed to make the results reproducible
set.seed(1)

covariates = all.vars(m$XFormula)[2:4]
covariatenames = c("Closterovirus T1", "Betapartitivirus T1", "P. lanceolata latent virus T1")
taxanames = c("CLOS", "ENAM", "BETA", "PILV", "CAUL")
PDFnames = c("PRED_prev_by_genotype_vs_ClosteroT1.pdf", "PRED_prev_by_genotype_vs_BetapartitiT1.pdf", "PRED_prev_by_genotype_vs_PilvT1.pdf")

for(n in 1:3)
{
  pdf(PDFnames[n], width = 10, height = 10)
  plots = list()
  
  covariate = covariates[n]
  covariatename = covariatenames[n]
  
  predictions_per_virus_T1 = data.frame(matrix(nrow = 4, ncol = 40))
  colnames(predictions_per_virus_T1) = paste0(sort(rep(taxanames, 8)), "_", covariate, c("_abs_G1", "_pres_G1", "_abs_G2", "_pres_G2", "_abs_G3", "_pres_G3", "_abs_G4", "_pres_G4"))
  rownames(predictions_per_virus_T1) = c("lo", "me", "hi", "probability")
  
  # For each host plant genotype (1-4) in turn:
  for(g in 1:4)
  {
    #And each virus taxon in turn:
    for(virus in 1:5)
    { 
      Gradient = constructGradient(m, focalVariable = covariate, non.focalVariables = 1)
      for(xlist in 1:5)
      {Gradient$XDataNew[[xlist]]$Genotype = factor(g, levels = c(1:4))}
      #Generate 1000 posterior predicted Y matrices (Predicted p(occurrences) per virus taxon):
      predY = predict(m, Gradient = Gradient, expected = TRUE)
      tmp = abind(predY, along = 3)
      #Extract the quantiles of the posterior predictions
      qpredY = apply(tmp, c(1, 2), quantile, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
      
      #Collate the predictions by focal virus, predictor virus T1 covariate and plant genotype
      taxon = taxanames[virus]
      
      lo = qpredY[1, , virus]
      predictions_per_virus_T1[1, paste(taxon, covariate, "abs", paste0("G", g), sep = "_")] = lo[1]
      predictions_per_virus_T1[1, paste(taxon, covariate, "pres", paste0("G", g), sep = "_")] = lo[2]
      hi = qpredY[3, , virus]
      predictions_per_virus_T1[3, paste(taxon, covariate, "abs", paste0("G", g), sep = "_")] = hi[1]
      predictions_per_virus_T1[3, paste(taxon, covariate, "pres", paste0("G", g), sep = "_")] = hi[2]
      me = qpredY[2, , virus]
      predictions_per_virus_T1[2, paste(taxon, covariate, "abs", paste0("G", g), sep = "_")] = me[1]
      predictions_per_virus_T1[2, paste(taxon, covariate, "pres", paste0("G", g), sep = "_")] = me[2]
      
      Pr = mean(tmp[2,virus,]>tmp[1,virus,])
      
      predictions_per_virus_T1[4, paste(taxon, covariate, "abs", paste0("G", g), sep = "_")] = Pr
      predictions_per_virus_T1[4, paste(taxon, covariate, "pres", paste0("G", g), sep = "_")] = Pr
      
      # Plot the predictions:
      linetype = ifelse(Pr>0.95|Pr<0.05, 1, 2)
      
      xx = Gradient$XDataNew[[virus]][,1]
      xlabel = colnames(Gradient$XDataNew[[1]])[1]
      ylabel = m$spNames[[virus]]
      main = paste0("Genotype ", g)
      
      toPlot = data.frame(xx, me, lo, hi, stringsAsFactors = TRUE)
      
      pl = ggplot(toPlot, aes(x = xx, y = me, fill = factor(xx))) + geom_bar(position = position_dodge(),
                                                                                    stat = "identity") + ylim(0, 1) + labs(title=main) + xlab(xlabel) + ylab(ylabel) +
        geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2,
                      position = position_dodge(0.9)) +
        theme(legend.position="none")
      
      plots[[length(plots) + 1]] = pl
    }
  }
  #Print all the plots into the PDF file named above:
  grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]], plots[[11]], plots[[12]], plots[[13]], plots[[14]], plots[[15]], plots[[16]], plots[[17]], plots[[18]], plots[[19]], plots[[20]], nrow = 4)
  dev.off()
  
  #Tidy up the predictions table and add columns defining the response variable, focal virus T1 covariate and genotype
  predictions_per_virus_T1 = data.frame(t(predictions_per_virus_T1))
  predictions_per_virus_T1$response = NA
  predictions_per_virus_T1$virusT1 = NA
  predictions_per_virus_T1$genotype = NA
  
  for(i in 1:nrow(predictions_per_virus_T1))
  {
    split_rowname = strsplit(rownames(predictions_per_virus_T1)[i], "_")[[1]]
    predictions_per_virus_T1$response[i] = split_rowname[1]
    predictions_per_virus_T1$virusT1[i] = paste(split_rowname[2], split_rowname[4], sep = "_")
    predictions_per_virus_T1$genotype[i] = split_rowname[5]
  }
  #Reorder data columns and export data as a csv file
  predictions_per_virus_T1 = predictions_per_virus_T1[,c(5:7, 1:4)]
  write.csv2(predictions_per_virus_T1, paste0("BarplotData_", covariates[n], ".csv"), row.names = FALSE)
}


          
