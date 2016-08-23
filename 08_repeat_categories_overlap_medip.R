# File: 08_repeat_categories_overlap_medip.R
# Auth: uhkniazi
# DESC: repeat categories that overlap with medip data
# Date: 23/08/2016


# source header file
source('00_header.R')
# jeffery's prior for gamma distribution
jef.prior = c(alpha=0.5, beta=0.0001)

### data loading
# connect to mysql database to find file locations
library('RMySQL')

##### connect to mysql database to get samples
db = dbConnect(MySQL(), user='rstudio', password='12345', dbname='Projects', host='127.0.0.1')
dbListTables(db)
dbListFields(db, 'MetaFile')
# another way to get the query, preferred
dfSample = dbGetQuery(db, "select * from MetaFile where idData=3;")
# close connection after getting data
dbDisconnect(db)

## pooled lifecycle stages with fisher p.values
loc = paste0(dfSample[dfSample$id == 2, 'location'], dfSample[dfSample$id == 2, 'name'])
oGRLpooled = f_LoadObject(loc)

loc = paste0(dfSample[dfSample$id == 1, 'location'], dfSample[dfSample$id == 1, 'name'])
oGRrep = f_LoadObject(loc)
## choose the appropriate feature for overlaps
oGRLCas = split(oGRrep, factor(oGRrep$fCategories))
head(names(oGRLCas))
length(oGRLCas)

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
# i = which(mCont['True',] > 0)
# mCont = mCont[,i]
temp = f_step2(mCont, main=paste(cvSample, '- proportion of Features with a MeDIP Signal'))
somule = temp

# female
cvSample = 'Female'
mCont = f_step1(oGRLCas, oGRLpooled$f)
temp = f_step2(mCont, main=paste(cvSample, '- proportion of Features with a MeDIP Signal'))
female = temp

# male
cvSample = 'Male'
mCont = f_step1(oGRLCas, oGRLpooled$m)
temp = f_step2(mCont, main=paste(cvSample, '- proportion of Features with a MeDIP Signal'))
male = temp

## perform a pairwise comparison for each category
lComparison = lapply(names(oGRLCas), function(x){
  df = data.frame(cbind(somule[,x], female[,x], male[,x]))
  st = stack(df)
  return(get.pairwise.prob(st$values, st$ind, p.adjust.method = 'bonf'))
})

ind = sapply(lComparison, function(x){
  any(na.omit(x <= 0.05))
})

lComparison.sig = lComparison[ind]
names(lComparison.sig) = names(oGRLCas)[ind]

## write the names of categories showing significant changes
dfResults = data.frame(category = names(lComparison.sig))

## plot the results for these categories
f_make.barplots = function(mDat, ...){
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
}

pdf('Results/significant.repeat.categories.pdf')
temp = sapply(names(lComparison.sig), function(x){
  df = data.frame(cbind(somule[,x], female[,x], male[,x]))
  colnames(df) = c('somule', 'female', 'male')
  f_make.barplots(df, main=x)
})
dev.off(dev.cur())

write.csv(dfResults, file='Results/significant.repeat.categories.csv')