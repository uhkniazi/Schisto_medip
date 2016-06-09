# File: 07_ncSeqs_overlaps_intragenic.R
# Auth: u.niazi@imperial.ac.uk
# DESC: uses the intragenic features overlapping with pooled medip data, repeats, ncSeqs and finds overlaps.
# Date: 09/06/2016


# source header file
source('00_header.R')
# jeffery's prior for gamma distribution
jef.prior = c(alpha=0.5, beta=0.0001)

### data loading
## pooled lifecycle stages with fisher p.values
oGRLpooled = f_LoadObject(file.choose())

## features object with ncSeqs
## load the repeats object
lFeatures = f_LoadObject(file.choose())
# sanity check
names(lFeatures)
sapply(lFeatures, length)
lFeatures$desc
lFeatures$desc.2
lReps = list(repeats=reduce(lFeatures$repeats), aber.miRNA=reduce(lFeatures$aber.miRNA),
             wormbase.rib=reduce(lFeatures$wormbase.rib), wormbase.nonRib=reduce(lFeatures$wormbase.nonRib)
             , mirBase=reduce(lFeatures$mirBase))

lCas = list(gene=reduce(lFeatures$gene), exon=reduce(lFeatures$exon), intron=reduce(lFeatures$intron),
            prom=reduce(lFeatures$prom), downstream=reduce(lFeatures$downstream))

oGRLCas = GRangesList(lCas)
names(lCas); sapply(lCas, length)
names(oGRLCas)

### steps for analysis
# step 1 - feature overlaps with a peak T/F
# step 2 - keep features if T in step 1
f_step1.2 = function(g, p){
  f = overlapsAny(g, p)
  return(g[f])
}

# step 3 - for each intragenic feature, find overlapping ncSeqs get 
#         contingency table of overlap
f_step3 = function(g, p){
  ## get the overlaps matrix
  lCounts = sapply(g, overlapsAny, p)
  lCounts = sapply(lCounts, function(x) factor(x, levels = c('TRUE', 'FALSE')))
  mDat = sapply(lCounts, table)
  rownames(mDat) = c('True', 'False')
  return(mDat)
}

## step 4 - plot the results in a bar plot
# calculate P(Feature | TRUE) i.e. distribution of the peaks over features as dirichlet distribution
f_step4 = function(mDat, ...){
  ## get dirichlet posterior with jeffery's prior
  # use a jeffery's dirichlet prior to calculate posterior
  alpha = mDat['True',] + 1/2
  if(!require(LearnBayes)) stop('Package LearnBayes required')
  mDat = rdirichlet(10000, alpha)
  colnames(mDat) = names(alpha)
  
  # get the median to plot
  mBar = apply(mDat, 2, mean)
  names(mBar) = colnames(mDat)
  l = barplot(mBar, beside=T, ylim=c(0,1), ...)
  ## draw error bars
  f_barplot_errorbars = function(x.loc, y.loc, ...){
    segments(x.loc, y.loc[1], x.loc, y.loc[2], ...)
    segments(x.loc-0.1, y.loc[1], x.loc+0.1, y.loc[1], ...)
    segments(x.loc-0.1, y.loc[2], x.loc+0.1, y.loc[2], ...)
  }
  sapply(seq_along(1:ncol(mDat)), function(x) f_barplot_errorbars(l[x,1], quantile(mDat[,x], c(0.025, 0.975))))
  return(mDat)
}

## step 5 - plot proportions 
# calculate P(TRUE | Feature) i.e. Proportion of features with a peak using beta distribution
## use a jeffery's prior for the proprtion parameter theta
getBetaPost = function(s, f){
  # simulate the theta for each data using a jeffery's prior beta prior
  rs = rbeta(10000, s+0.5, f+0.5)
  return(rs)  
}

f_step5 = function(mDat, ...){
  ## get the beta posterior across columns
  cn = colnames(mDat)
  suc = mDat['True',]
  fail = mDat['False',]
  # get the posterior distribution of the parameter theta
  mDat = sapply(seq_along(1:length(suc)), function(x){
    return(getBetaPost(suc[x], fail[x]))
  })
  
  colnames(mDat) = cn
  # get the mean to plot
  mBar = apply(mDat, 2, mean)
  names(mBar) = colnames(mDat)
  m = apply(mDat, 2, function(x) quantile(x, 0.98))
  l = barplot(mBar, beside=T, ylim=c(0,max(m)), ...)
  ## draw error bars
  f_barplot_errorbars = function(x.loc, y.loc, ...){
    segments(x.loc, y.loc[1], x.loc, y.loc[2], ...)
    segments(x.loc-0.1, y.loc[1], x.loc+0.1, y.loc[1], ...)
    segments(x.loc-0.1, y.loc[2], x.loc+0.1, y.loc[2], ...)
  }
  sapply(seq_along(1:ncol(mDat)), function(x) f_barplot_errorbars(l[x,1], quantile(mDat[,x], c(0.025, 0.975))))
  return(mDat)
}

f_getSummary = function(x){
  return(c(mean=mean(x), std.dev=sd(x), median=median(x), ci.low=quantile(x, 0.025), ci.high=quantile(x, 0.975)))
}


## repeat the process for each lifecycle stage
# somules
lOverlap = lapply(oGRLCas, f_step1.2, oGRLpooled$s)
cvNames = names(lOverlap)
cvSample = 'Somule'
temp = sapply(seq_along(cvNames), function(x){
  cTitle = cvNames[x]
  mCont = f_step3(lReps, lOverlap[[cTitle]])
  f_step4(mCont, main=paste(cvSample, '- distribution of peaks over ncSeqs overlapping with', cTitle))
  f_step5(mCont, main=paste(cvSample, '- proportion of ncSeqs with peak and overlapping with', cTitle))
  print(cTitle)
  print(mCont)
})


lOverlap = lapply(oGRLCas, f_step1.2, oGRLpooled$f)
cvNames = names(lOverlap)
cvSample = 'Female'
temp = sapply(seq_along(cvNames), function(x){
  cTitle = cvNames[x]
  mCont = f_step3(lReps, lOverlap[[cTitle]])
  f_step4(mCont, main=paste(cvSample, '- distribution of peaks over ncSeqs overlapping with', cTitle))
  f_step5(mCont, main=paste(cvSample, '- proportion of ncSeqs with peak and overlapping with', cTitle))
  print(cTitle)
  print(mCont)
})


lOverlap = lapply(oGRLCas, f_step1.2, oGRLpooled$m)
cvNames = names(lOverlap)
cvSample = 'Male'
temp = sapply(seq_along(cvNames), function(x){
  cTitle = cvNames[x]
  mCont = f_step3(lReps, lOverlap[[cTitle]])
  f_step4(mCont, main=paste(cvSample, '- distribution of peaks over ncSeqs overlapping with', cTitle))
  f_step5(mCont, main=paste(cvSample, '- proportion of ncSeqs with peak and overlapping with', cTitle))
  print(cTitle)
  print(mCont)
})

