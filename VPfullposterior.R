#.
#.              FULL POSTERIOR DISTRIBUTION - VARIANCE PARTITIONING
#.              ---------------------------------------------------
#.
#... The below code extends the variance partitioning R function available in github at https://github.com/hmsc-r/HMSC
#... called computeVariancePartitioning from the HMSC package (Tikhonov et al., 2022). Downloaded in June 2022.
#... This code focuses only on variance partition itself, thus other computations included in the original function have been 
#... disregarded and commented out of the script.
#... 
#... Code extended by Emy Guilbault, September 2022. Last changes June 2023.
#... 
#... 
#..................................................................................

library(Hmsc)

load('models_thin_5000_samples_250_chains_4.Rdata')
# name of the object of interest: models


### modified function to calculate the full VP posterior distribution, rather than mean VP fractions only:
computeVariancePartitioningfull = function(hM, group=NULL, groupnames=NULL, start=1, na.ignore=FALSE)
{
  
  # Name and prepare parameters 
  ny = hM$ny
  nc = hM$nc
  nt = hM$nt
  ns = hM$ns
  np = hM$np
  nr = hM$nr
  
  if(is.null(group)){
    ## default: use terms
    if(nc > 1){
      if (is.null(hM$XFormula))
        stop("no XFormula: you must give 'group' and 'groupnames'")
      group = attr(hM$X, "assign")
      
      if (class(hM$X) == "list")
        group = attr(hM$X[[1]], "assign")
        
      if (group[1] == 0) # assign (Intercept) to group
        group[1] <- 1
      groupnames = attr(terms(hM$XFormula, data = hM$XData), "term.labels")
    } else {
      group = c(1)
      groupnames = hM$covNames[1]
    }
    
  }
  
  
  ngroups = max(group)
  
  
  R2T.Y = 0
  R2T.Beta = rep(0,nc)
  
  
  ## Retrieve info over MCMC iterations
  postList=poolMcmcChains(hM$postList, start=start)
  
  
  # +++++++ prepare to save all info in an array 
  # Initialization of the arrays
  nsamp = length(postList) # to get number of samples for all chains (this should be now corrected in the guthub code as well)
  
  fixed = array(0, dim=c(ns, 1, nsamp))
  fixedsplit = array(0, dim=c(ns, ngroups, nsamp))
  random = array(0, c(ns, nr, nsamp))
  
  
  
  #If na.ignore=T, convert XData to a list
  if(na.ignore){
    xl=list()
    for(s in 1:ns){
      xl[[s]]=hM$X
    }
    hM$X=xl
  }
  
  switch(class(hM$X)[1L], matrix = {
    cMA = cov(hM$X)
  }, list = {
    
    if(na.ignore){
      cMA = list()
      for(s in 1:ns){cMA[[s]]=cov(hM$X[[s]][which(hM$Y[,s]>-Inf),])}
    }
    else{cMA = lapply(hM$X, cov)}
  })
  
  
  #***** The following part of the original HMSC computeVariancePartitioning code is not used here, in calculating the full posterior distribution of the VP
  # it is thus disregarded here. Similarly the lines commented out are tagged with ***
  
  #*** geta = function(a){
  # switch(class(hM$X)[1L],
  # matrix = {res = hM$X%*%(t(hM$Tr%*%t(a$Gamma)))},
  # list = {
  # res = matrix(NA,hM$ny,hM$ns)
  # for(j in 1:hM$ns) res[,j] =  hM$X[[j]]%*%(t(hM$Tr[j,]%*%t(postList[[1]]$Gamma)))
  # }
  # )
  # return(res)
  # }
  # la=lapply(postList, geta)
  
  # getf = function(a){
  # switch(class(hM$X)[1L],
  # matrix = {res = hM$X%*%(a$Beta)},
  # list = {
  # res = matrix(NA,hM$ny,hM$ns)
  # for(j in 1:hM$ns) res[,j] = hM$X[[j]]%*%a$Beta[,j]
  # }
  # )
  # return(res)
  # }
  # lf=lapply(postList, getf)
  
  # gemu = function(a){
  # res = t(hM$Tr%*%t(a$Gamma))
  # return(res)
  # }
  # lmu=lapply(postList, gemu)
  
  
  # gebeta = function(a){
  # res = a$Beta
  # return(res)
  # }
  # lbeta=lapply(postList, gebeta
  
  
  
  ##  run over MCMC chains samples 
  for (i in 1:nsamp){
    #***
    #for (k in 1:nc){
    #   R2T.Beta[k] = R2T.Beta[k] + cor(lbeta[[i]][k,],lmu[[i]][k,])^2
    # }
    
    fixedsplit1 = matrix(0 ,nrow=ns, ncol=ngroups)
    fixed1 = matrix(0, nrow=ns, ncol=1)
    
    random1 = matrix(0, nrow=ns, ncol=nr) 
    Beta = postList[[i]]$Beta
    Lambdas = postList[[i]]$Lambda
    
    #***
    #f = lf[[i]]
    #a = la[[i]]
    
    #a=a-matrix(rep(rowMeans(a),hM$ns),ncol=hM$ns)
    #f=f-matrix(rep(rowMeans(f),hM$ns),ncol=hM$ns)
    #res1 = sum((rowSums((a*f))/(hM$ns-1))^2)
    #res2 = sum((rowSums((a*a))/(hM$ns-1))*(rowSums((f*f))/(hM$ns-1)))
    #R2T.Y = R2T.Y + res1/res2
    for (j in 1:ns){
      switch(class(hM$X)[1L], matrix = {cM = cMA}, list = {cM = cMA[[j]]})
      
      ftotal = t(Beta[,j])%*%cM%*%Beta[,j]
      fixed1[j] = fixed1[j] + ftotal
      for (k in 1:ngroups){
        sel = (group==k)
        fpart = t(Beta[sel,j])%*%cM[sel,sel]%*%Beta[sel,j]
        fixedsplit1[j,k] = fixedsplit1[j,k] + fpart
      }
      
    }
    
    for (level in seq_len(nr)){
      Lambda = Lambdas[[level]]
      nf = dim(Lambda)[[1]]
      for (factor in 1:nf){
        random1[,level] = random1[,level] + t(Lambda[factor,])*Lambda[factor,]
      }
    }
    
    if (nr>0){
      tot = fixed1 + rowSums(random1)
      fixed[,,i] = fixed[,,i] + fixed1/tot
      for (level in seq_len(nr)){
        random[,level,i] = random[,level,i] + random1[,level]/tot
      }
    }else{
      fixed[,,i] = fixed[,,i] + matrix(1,nrow=ns,ncol=1)
    }
    
    for (k in 1:ngroups){
      fixedsplit[,k,i] =  fixedsplit[,k,i] + fixedsplit1[,k]/rowSums(fixedsplit1)
    }
    
    
  }
  
  fixed = fixed#/nsamp
  random = random#/nsamp
  fixedsplit = fixedsplit#/nsamp
  
  vals = array(0, dim=c(ngroups+nr, ns, nsamp))
  for (i in 1:nsamp) {
    for (j in 1:ngroups){
      vals[j,,i] = fixed[,,i]*fixedsplit[,j,i]
      
    }
    for (j in seq_len(nr)){
      vals[ngroups+j,,i] = random[,j,i]
    }
  }
  
  #***
  #R2T.Y = R2T.Y/hM$samples
  #R2T.Beta = R2T.Beta/hM$samples
  
  
  names(R2T.Beta)=hM$covNames
  colnames(vals)=hM$spNames
  leg = groupnames
  for (r in seq_len(nr)) {
    leg = c(leg, paste("Random: ", hM$rLNames[r], sep = ""))
  }
  
  rownames(vals)=leg

  
  VP = list()
  VP$vals = vals
  #***
  #VP$R2T = list(Beta=R2T.Beta,Y=R2T.Y)
  VP$group = group
  VP$groupnames = groupnames
  
  return(VP)
}

VPfull = computeVariancePartitioningfull(hM=models[[1]])
dim(VPfull$vals) #  Covarariates / Species  / MCMC samples

Mean = round(apply(VPfull$vals, c(1,2), mean), 3)
colSums(Mean)

# calculate a credible interval
CI.low = round(apply(VPfull$vals, c(1,2), function(x){quantile(x, 0.025)}), 3)
CI.up = round(apply(VPfull$vals, c(1,2),function(x){quantile(x, 0.975)}), 3)

colSums(CI.low)
colSums(CI.up)
