library(Hmsc)
library(ggplot2)


load("MJ_models_thin_5000_samples_250_chains_4.Rdata")
m = models[[1]]
mpost = convertToCodaObject(m)
CIs = summary(mpost$Beta)
betatable = data.frame(rownames(CIs[[1]]),CIs[[1]][, 1], CIs[[2]][, c(3, 1, 5)])
names(betatable) = c("parameter","mean", "median", "lower_CI", "upper_CI")
betatable$dataset <- "pooled"


load("MJ_T2_models_thin_5000_samples_250_chains_4.Rdata")
m = models[[1]]
mpost = convertToCodaObject(m)
CIs = summary(mpost$Beta)
betatableT = data.frame(rownames(CIs[[1]]),CIs[[1]][, 1], CIs[[2]][, c(3, 1, 5)])
names(betatableT) = c("parameter","mean", "median", "lower_CI", "upper_CI")
betatableT$dataset <- "T2"
betatable <- rbind(betatable,betatableT)

load("MJ_T3_models_thin_5000_samples_250_chains_4.Rdata")
m = models[[1]]
mpost = convertToCodaObject(m)
CIs = summary(mpost$Beta)
betatableT = data.frame(rownames(CIs[[1]]),CIs[[1]][, 1], CIs[[2]][, c(3, 1, 5)])
names(betatableT) = c("parameter","mean", "median", "lower_CI", "upper_CI")
betatableT$dataset <- "T3"
betatable <- rbind(betatable,betatableT)

load("MJ_T4_models_thin_5000_samples_250_chains_4.Rdata")
m = models[[1]]
mpost = convertToCodaObject(m)
CIs = summary(mpost$Beta)
betatableT = data.frame(rownames(CIs[[1]]),CIs[[1]][, 1], CIs[[2]][, c(3, 1, 5)])
names(betatableT) = c("parameter","mean", "median", "lower_CI", "upper_CI")
betatableT$dataset <- "T4"
betatable <- rbind(betatable,betatableT)

parametersplit <- unlist(strsplit(betatable$parameter," "))

betatable$effect <- substr(parametersplit[seq(1,length(parametersplit),by=4)],start = 3,stop = 20)
betatable$response <- parametersplit[seq(3,length(parametersplit),by=4)]

saveRDS(betatable,"betatable.rds")

betatable$response1<-factor(betatable$response, levels=c("CLOSTERO", "BETAPARTITI", "PILV", "CAULIMO", "ENAMO"))

pdf("beta_interval.pdf", width = 14, height = 9)
p <- ggplot(betatable,aes(x=median,xmin=lower_CI, xmax = upper_CI,y=effect,color=dataset,shape=dataset))+
  geom_pointrange(linewidth=0.6,position = position_dodge(width = 0.6)) + 
  facet_grid(cols=vars(response1)) + geom_vline(xintercept=0,linetype="dotted") +
  xlim(-10,5)
print(p2)
p

dev.off()