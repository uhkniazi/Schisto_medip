# File: 07_features_overlap_medip.R
# Auth: u.niazi@imperial.ac.uk
# DESC: generic script to calculate overlaps of features (query), and medip (subject)
# Date: 01/07/2016


# source header file
source('00_header.R')
# jeffery's prior for gamma distribution
jef.prior = c(alpha=0.5, beta=0.0001)

### data loading
## pooled lifecycle stages with fisher p.values
oGRLpooled = f_LoadObject(file.choose())

## features object lFeatures
lFeatures = f_LoadObject(file.choose())
# sanity check
names(lFeatures)
sapply(lFeatures, length)
lFeatures$desc.3
lReps = list(repeats=reduce(lFeatures$repeats), mirbase=reduce(lFeatures$mirbase),
             wormbase=reduce(lFeatures$wormbase))

lCas = list(gene=reduce(lFeatures$gene), exon=reduce(lFeatures$exon), intron=reduce(lFeatures$intron),
            prom=reduce(lFeatures$prom), downstream=reduce(lFeatures$downstream))

lGenomic = list(intra=reduce(lFeatures$intra), inter=reduce(lFeatures$inter))

## choose the appropriate feature for overlaps
oGRLCas = GRangesList(lGenomic)
sapply(oGRLCas, length)
names(oGRLCas)

### steps for analysis
# step 1 - find overlaps of feature (query) with medip (subject)
f_step1 = function(q, s){
  ## get the overlaps matrix
  lCounts = sapply(q, overlapsAny, s)
  lCounts = sapply(lCounts, function(x) factor(x, levels = c('TRUE', 'FALSE')))
  mDat = sapply(lCounts, table)
  rownames(mDat) = c('True', 'False')
  return(mDat)
}

# ## step 4 - plot the results in a bar plot
# # calculate P(Feature | TRUE) i.e. distribution of the peaks over features as dirichlet distribution
# f_step4 = function(mDat, ...){
#   ## get dirichlet posterior with jeffery's prior
#   # use a jeffery's dirichlet prior to calculate posterior
#   alpha = mDat['True',] + 1/2
#   if(!require(LearnBayes)) stop('Package LearnBayes required')
#   mDat = rdirichlet(10000, alpha)
#   colnames(mDat) = names(alpha)
#   
#   # get the median to plot
#   mBar = apply(mDat, 2, mean)
#   names(mBar) = colnames(mDat)
#   l = barplot(mBar, beside=T, ylim=c(0,1), ...)
#   ## draw error bars
#   f_barplot_errorbars = function(x.loc, y.loc, ...){
#     segments(x.loc, y.loc[1], x.loc, y.loc[2], ...)
#     segments(x.loc-0.1, y.loc[1], x.loc+0.1, y.loc[1], ...)
#     segments(x.loc-0.1, y.loc[2], x.loc+0.1, y.loc[2], ...)
#   }
#   sapply(seq_along(1:ncol(mDat)), function(x) f_barplot_errorbars(l[x,1], quantile(mDat[,x], c(0.025, 0.975))))
#   return(mDat)
# }

## step 2 - plot proportions 
# calculate P(TRUE | Feature) i.e. Proportion of features with a peak using beta distribution
## use a jeffery's prior for the proprtion parameter theta
getBetaPost = function(s, f){
  # simulate the theta for each data using a jeffery's prior beta prior
  rs = rbeta(10000, s+0.5, f+0.5)
  return(rs)  
}

f_step2 = function(mDat, ...){
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
cvSample = 'Somule'
mCont = f_step1(oGRLCas, oGRLpooled$s)
temp = f_step2(mCont, main=paste(cvSample, '- proportion of Features with a MeDIP Signal'))
print(cvSample)
print(mCont)
apply(temp, 2, function(x) format(f_getSummary(x), digits = 3))

# female
cvSample = 'Female'
mCont = f_step1(oGRLCas, oGRLpooled$f)
temp = f_step2(mCont, main=paste(cvSample, '- proportion of Features with a MeDIP Signal'))
print(cvSample)
print(mCont)
apply(temp, 2, function(x) format(f_getSummary(x), digits = 3))

# male
cvSample = 'Male'
mCont = f_step1(oGRLCas, oGRLpooled$m)
temp = f_step2(mCont, main=paste(cvSample, '- proportion of Features with a MeDIP Signal'))
print(cvSample)
print(mCont)
apply(temp, 2, function(x) format(f_getSummary(x), digits = 3))


