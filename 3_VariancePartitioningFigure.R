setwd()

load("models/models_thin_5000_samples_250_chains_4.Rdata")
m = models[[1]]

library(Hmsc)

pdf("VarPart.pdf", width = 6, height = 5)
par(mar = c(5, 5, 3, 8.5), lwd = 0.1)
{
  preds = computePredictedValues(m)
  VP = computeVariancePartitioning(m, group = attr(m$X[[1]], "assign"), groupnames = attr(terms(m$XFormula), "term.labels"))
  mycols = rainbow(nrow(VP$vals))
  MF = evaluateModelFit(hM=m, predY=preds)
  R2 = MF$TjurR2
  barplot(VP$vals, ylim = c(0,1), cex.lab = 0.8, ylab = "Fractions of variance explained", cex.axis = 0.7, col = mycols, las=2, horiz=FALSE, cex.names = 0.6)
  legend("topright", inset = c(-0.48, 0), xpd=TRUE, legend = rev(c(names(m$XData[[1]]), "Random (plant)", "Random (population)")), fill = rev(mycols), bty = "n", cex = 0.8)
}
dev.off()