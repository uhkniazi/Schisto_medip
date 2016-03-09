# File: 02_features_overlaps_casette_3.R
# Auth: u.niazi@imperial.ac.uk
# DESC: uses the casette object and pooled medip peaks for each lifecycle stage created earlier, 
#       creates a gene casette and finds overlaps of casettes and medip peaks. Reports statistics and plots
# Date: 01/05/2015


# source header file
source('00_header.R')
# jeffery's prior for gamma distribution
jef.prior = c(alpha=0.5, beta=0.0001)

### data loading
## pooled lifecycle stages with fisher p.values
oGRLpooled = f_LoadObject(file.choose())

## features object with gene casette information created earlier
lFeatures = f_LoadObject(file.choose())

# sanity check
sapply(lFeatures, length)
lFeatures$desc
sapply(oGRLpooled, length)
metadata(oGRLpooled)


lCas = lFeatures
lCas$desc = NULL

### steps for analysis
# step 1 - peak intersects with a gene , TRUE or FALSE
# step 2 - use peaks TRUE at step 1
# step 3 - if it intersects with more than one gene then discard
# step 4 - assign gene id to peak. 
# step 5 - get contingency table
# step 6 - calculate P(Feature | TRUE) i.e. distribution of the peaks over features as dirichlet distribution
# step 7 - calculate P(TRUE | Feature) i.e. Proportion of features with a peak using beta distribution

### functions used in script
## use a jeffery's prior for the proprtion parameter theta
getBetaPost = function(s, f){
  # simulate the theta for each data using a jeffery's prior beta prior
  rs = rbeta(10000, s+0.5, f+0.5)
  return(rs)  
}

f_step1 = function(g, p){
  # the range over which to find the overlaps
  g1 = unlist(g$gene)
  g2 = unlist(g$upstream)
  g3 = unlist(g$downstream)
  g.all = reduce(GRangesList(g1, g2, g3))
  f = overlapsAny(p, g.all)
  return(f)
}

f_step3 = function(g, p){
  # check if a peak overlaps more than one genes
  c = countOverlaps(p, g$gene)
  i = which(c > 1)
  p = p[-i]
  # check if peak overlaps more than one upstream or downstream regions
  c = countOverlaps(p, g$upstream)
  i = which(c > 1)
  p = p[-i]
  c = countOverlaps(p, g$downstream)
  i = which(c > 1)
  p = p[-i]
  return(p)
}


f_step4 = function(g, p){
  p$GeneID = ''
  df = as.data.frame(findOverlaps(p, g$gene))
  n = names(g$gene)
  p[df$queryHits]$GeneID = n[df$subjectHits]
  
  # upstream region
  df = as.data.frame(findOverlaps(p, g$upstream))
  n = names(g$upstream)
  p$UpStream = ''
  # some peaks can intersect a gene and upstream region
  p[df$queryHits]$UpStream = n[df$subjectHits]
  
  # downstream
  df = as.data.frame(findOverlaps(p, g$downstream))
  n = names(g$downstream)
  p$DownStream = ''
  # some peaks can intersect a gene and upstream region
  p[df$queryHits]$DownStream = n[df$subjectHits]
  return(p)
}

f_step5 = function(g, p){
  ## get the overlaps matrix
    us = overlapsAny(unlist(g$upstream), p)
    ds = overlapsAny(unlist(g$downstream), p)
    fst.exon = overlapsAny(unlist(g$fst.exon), p)
    exons.ot = overlapsAny(unlist(g$exons.others), p)
    introns = overlapsAny(unlist(g$introns), p)
    
    mDat = cbind(us=as.numeric(table(us)), fst.exon=as.numeric(table(fst.exon)), exons.ot=as.numeric(table(exons.ot)),
                 introns=as.numeric(table(introns)), ds=as.numeric(table(ds)))
    
    rownames(mDat) = c('False', 'True')
    return(mDat)
}


f_step6 = function(mDat, ...){
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


f_step7 = function(mDat, ...){
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
  l = barplot(mBar, beside=T, ylim=c(0,max(mDat)), ...)
  ## draw error bars
  f_barplot_errorbars = function(x.loc, y.loc, ...){
    segments(x.loc, y.loc[1], x.loc, y.loc[2], ...)
    segments(x.loc-0.1, y.loc[1], x.loc+0.1, y.loc[1], ...)
    segments(x.loc-0.1, y.loc[2], x.loc+0.1, y.loc[2], ...)
  }
  sapply(seq_along(1:ncol(mDat)), function(x) f_barplot_errorbars(l[x,1], quantile(mDat[,x], c(0.025, 0.975))))
  return(mDat)
}

# f_step5 = function(g, p){
#   # create a matrix to hold the matches
#   mFeatures = matrix(NA, nrow = length(p), ncol = 5)
#   colnames(mFeatures) = c('upstream', 'fst.exon', 'exons.others', 'introns', 'downstream')
#   cn = colnames(mFeatures)
#   # fill the data
#   for(i in 1:ncol(mFeatures)){
#     # fill data one column at a time
#     mFeatures[,i] = countOverlaps(p, g[[cn[i]]])
#   }
#   df = as.data.frame(mcols(p))
#   df = cbind(df, mFeatures)
#   mcols(p) = df
#   return(list(peaks=p, overlaps=mFeatures))
# }
# 
# f_step6 = function(s5, ...){
#   ## summarizing the proprtion of the presense or absence of peaks
#   suc = colSums((s5$overlaps))
#   total = nrow(s5$overlaps)
#   fail = total - suc
#   # get the posterior distribution of the parameter theta
#   mDat = sapply(seq_along(1:length(suc)), function(x){
#     return(getBetaPost(suc[x], fail[x]))
#   })
#   
#   colnames(mDat) = colnames(s5$overlaps)
#   # get the mean to plot
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

# f_step7 = function(g, p, ...){
# #   boot.rate = function(qu, sub){
# #     iBootnum = 100
# #     mBoot = matrix(NA, iBootnum, 2)
# #     colnames(mBoot) = c('False', 'True')
# #     for (i in 1:iBootnum){
# #       sam = sample(1:length(qu), size = 1000, replace = T)
# #       mBoot[i,] = as.numeric(table(factor(overlapsAny(qu[sam], sub), levels = c('FALSE', 'TRUE'))))
# #     }
# #     return((mBoot))
# #   }
# # uncomment this section as not using bootstrap to calculate rate
#   us = overlapsAny(unlist(g$upstream), p)
#   ds = overlapsAny(unlist(g$downstream), p)
#   fst.exon = overlapsAny(unlist(g$fst.exon), p)
#   exons.ot = overlapsAny(unlist(g$exons.others), p)
#   introns = overlapsAny(unlist(g$introns), p)
#   
#   mDat = cbind(us=as.numeric(table(us)), fst.exon=as.numeric(table(fst.exon)), exons.ot=as.numeric(table(exons.ot)),
#                introns=as.numeric(table(introns)), ds=as.numeric(table(ds)))
#   
#   rownames(mDat) = c('False', 'True')
#   ivDat = colSums(mDat)
#   # convert to rates
#   ivProb = sweep(mDat, 2, ivDat, '/')['True',] * 1000
#   ## if not using bootstrap then comment this section and uncomment previous
#   us = boot.rate(unlist(g$upstream), p)
#   ds = boot.rate(unlist(g$downstream), p)
#   fst.exon = boot.rate(unlist(g$fst.exon), p)
#   exons.ot = boot.rate(unlist(g$exons.others), p)
#   introns = boot.rate(unlist(g$introns), p)
#   
#   mDat = cbind(us=us[,'True'], fst.exon=fst.exon[,'True'], exons.ot=exons.ot[,'True'],
#                introns=introns[,'True'], ds=ds[,'True'])
#   #rownames(mDat) = c('False', 'True')
#   # convert to rates
#   #ivProb = mDat['True',]
#   
#   #prior = rbeta(1000, ivDat['True'], ivDat['False'])
#   mDat = apply(mDat, 2, function(x) {
#     r = rgamma(1000, shape = jef.prior['alpha'] + sum(x) , jef.prior['beta']+ length(x))
#     #r = rgamma(1000, shape = jef.prior['alpha'] + (x*100) , jef.prior['beta']+ 100)
#     return(r)
#   })
#   
#   mBar = apply(mDat, 2, mean)
#   names(mBar) = colnames(mDat)
#   l = barplot(mBar, beside=T, ylim=c(0,max(mDat)), ...)
#   ## draw error bars
#   f_barplot_errorbars = function(x.loc, y.loc, ...){
#     segments(x.loc, y.loc[1], x.loc, y.loc[2], ...)
#     segments(x.loc-0.1, y.loc[1], x.loc+0.1, y.loc[1], ...)
#     segments(x.loc-0.1, y.loc[2], x.loc+0.1, y.loc[2], ...)
#   }
#   sapply(seq_along(1:ncol(mDat)), function(x) f_barplot_errorbars(l[x,1], quantile(mDat[,x], c(0.025, 0.975))))
#   return(mDat)
# }
# 
# f_step8 = function(g, p, ...){
#   us = overlapsAny(unlist(g$upstream), p)
#   ds = overlapsAny(unlist(g$downstream), p)
#   fst.exon = overlapsAny(unlist(g$fst.exon), p)
#   exons.ot = overlapsAny(unlist(g$exons.others), p)
#   introns = overlapsAny(unlist(g$introns), p)
#   
#   mDat = cbind(us=as.numeric(table(us)), fst.exon=as.numeric(table(fst.exon)), exons.ot=as.numeric(table(exons.ot)),
#                introns=as.numeric(table(introns)), ds=as.numeric(table(ds)))
#   
#   rownames(mDat) = c('False', 'True')
#   
#   # use a jeffery's dirichlet prior to calculate posterior
#   alpha = mDat['True',] + 1/2
#   if(!require(LearnBayes)) stop('Package LearnBayes required')
#   mDat = rdirichlet(1000, alpha)
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

## load the repeats object
lFeatures = f_LoadObject(file.choose())
# sanity check
names(lFeatures)
sapply(lFeatures, length)
lFeatures$desc
### note: use this line if removing repeats
lCas = lapply(lCas, function(x) setdiff(unlist(x), lFeatures$repeats))
sapply(lCas, length)

# somules
f = f_step1(lCas, oGRLpooled$s)
s = f_step3(lCas, oGRLpooled$s[f])
s = f_step4(lCas, s)
mat = f_step5(lCas, s)
mDir.som = f_step6(mat, main='Somule - distribution of peaks over features')
mBeta.som = f_step7(mat, main='Somule - Proportion of features with a peak')
somule = s

summary(signif(mDir.som, 2))
summary(signif(mBeta.som, 2))

# female
f = f_step1(lCas, oGRLpooled$f)
s = f_step3(lCas, oGRLpooled$f[f])
s = f_step4(lCas, s)
mat = f_step5(lCas, s)
mDir.fem = f_step6(mat, main='Female - distribution of peaks over features')
mBeta.fem = f_step7(mat, main='Female - Proportion of features with a peak')
female = s

summary(signif(mDir.fem, 2))
summary(signif(mBeta.fem, 2))

# male
f = f_step1(lCas, oGRLpooled$m)
s = f_step3(lCas, oGRLpooled$m[f])
s = f_step4(lCas, s)
mat = f_step5(lCas, s)
mDir.male = f_step6(mat, main='Male - distribution of peaks over features')
mBeta.male = f_step7(mat, main='Male - Proportion of features with a peak')
male = s

summary(signif(mDir.male, 2))
summary(signif(mBeta.male, 2))


# plot them together
# get the mean to plot
## set the appropriate variables on what we want to plot
## just a quick hack to save from writing
mSom = mDir.som
mFem = mDir.fem
mMale = mDir.male
## OR
mSom = mBeta.som
mFem = mBeta.fem
mMale = mBeta.male

s = apply(mSom, 2, mean)
f = apply(mFem, 2, mean)
m = apply(mMale, 2, mean)

mBar = rbind(s, f, m)
col = grey.colors(nrow(mBar))
l = barplot(mBar, beside=T, ylim=c(0,max(mMale)), col=col, main='Proportion of features with a peak')
## draw error bars
f_barplot_errorbars = function(x.loc, y.loc, ...){
  segments(x.loc, y.loc[1], x.loc, y.loc[2], ...)
  segments(x.loc-0.1, y.loc[1], x.loc+0.1, y.loc[1], ...)
  segments(x.loc-0.1, y.loc[2], x.loc+0.1, y.loc[2], ...)
}
sapply(seq_along(1:ncol(mBar)), function(x) f_barplot_errorbars(l[1,x], quantile(mSom[,x], c(0.025, 0.975))))
sapply(seq_along(1:ncol(mBar)), function(x) f_barplot_errorbars(l[2,x], quantile(mFem[,x], c(0.025, 0.975))))
sapply(seq_along(1:ncol(mBar)), function(x) f_barplot_errorbars(l[3,x], quantile(mMale[,x], c(0.025, 0.975))))

legend('topleft', legend = c(rownames(mBar)), fill=col)


oGRLsamples.fisher.pooled = GRangesList(somule, female, male)
names(oGRLsamples.fisher.pooled) = c('somule', 'female', 'male')
sapply(oGRLsamples.fisher.pooled, length)

## create a gff for each sample
for (i in 1:length(oGRLsamples.fisher.pooled)){
  r = oGRLsamples.fisher.pooled[[i]]
  n = make.names(paste('peaks_overlapping_with_casette', names(oGRLsamples.fisher.pooled)[i], date(), 'gff'))
  dir.create('Results_scripts/Gff_files', showWarnings = F)
  export(r, con=paste('Results_scripts/Gff_files/', n, sep=''), format = 'gff3')
  n = make.names(paste('peaks_overlapping_with_casette', names(oGRLsamples.fisher.pooled)[i], date(), 'csv'))
  write.csv(as.data.frame(r), file=paste('Results_scripts/', n, sep=''))  
}

