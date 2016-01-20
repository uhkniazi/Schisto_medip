# File: 03_02_genes_ranking_over_conserved_peaks.R
# Auth: u.niazi@imperial.ac.uk
# DESC: uses the casette object and conserved medip peaks for across lifecycle stages created earlier, 
#       and ranks the genes based on overlapping pvalue
# Date: 19/01/2016


# source header file
source('00_header.R')

### data loading
## conserved lifecycle stages with fisher p.values
oGRconserved = f_LoadObject(file.choose())

## features object with gene casette information created earlier
lFeatures = f_LoadObject(file.choose())

# sanity check
sapply(lFeatures, length)
lFeatures$desc
metadata(oGRconserved)

lCas = lFeatures
lCas$desc = NULL

### steps for analysis
# step 1 - peak intersects with a gene , TRUE or FALSE
# step 2 - use peaks TRUE at step 1
# step 3 - if it intersects with more than one gene then discard
# step 4 - assign gene id to peak. 
# step 5 - which feature the peak intersects with, T/F for each feature
# step 6 - plot the posterior beta distribution for the binomial parameter theta

### functions used in script
## use a non-informative prior for the proprtion parameter theta
getPost = function(s, f){
  # simulate the theta for each data using a uniform beta prior
  rs = rbeta(1000, s+1, f+1)
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
  # create a matrix to hold the matches
  mFeatures = matrix(NA, nrow = length(p), ncol = 5)
  colnames(mFeatures) = c('upstream', 'fst.exon', 'exons.others', 'introns', 'downstream')
  cn = colnames(mFeatures)
  # fill the data
  for(i in 1:ncol(mFeatures)){
    # fill data one column at a time
    mFeatures[,i] = countOverlaps(p, g[[cn[i]]])
  }
  df = as.data.frame(mcols(p))
  df = cbind(df, mFeatures)
  mcols(p) = df
  return(list(peaks=p, overlaps=mFeatures))
}

f_step6 = function(s5, ...){
  ## summarizing the proprtion of the presense or absence of peaks
  suc = colSums((s5$overlaps))
  total = nrow(s5$overlaps)
  fail = total - suc
  # get the posterior distribution of the parameter theta
  mDat = sapply(seq_along(1:length(suc)), function(x){
    return(getPost(suc[x], fail[x]))
  })
  
  colnames(mDat) = colnames(s5$overlaps)
  # get the median to plot
  mBar = apply(mDat, 2, median)
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


f = f_step1(lCas, oGRconserved)
table(f)
s = f_step3(lCas, oGRconserved[f])
length(s)
s = f_step4(lCas, s)
head(s)
s = f_step5(lCas, s)
# plot the results
mOverlap = f_step6(s, main='Overlaps')
oGRpeaks = s$peaks
head(oGRpeaks)
# for each peak get the associated gene id as the splitting factor
fSplit = rep(NA, length.out=length(oGRpeaks))
# save unusual peaks over multiple genes
iUnusual = NA
for (i in seq_along(fSplit)){
  x = unique(c(oGRpeaks$GeneID[i], oGRpeaks$UpStream[i], oGRpeaks$DownStream[i]))
  fSplit[i] = x[x != ''][1]
  if (length(x[x!='']) > 1) iUnusual = append(iUnusual, i);
}

oGRLpeaks = split(oGRpeaks, fSplit)

oGRgenes.p = GRanges()

for (i in 1:length(oGRLpeaks)){
  gr = oGRLpeaks[[i]]
  pos = which.min(gr$fisher.p.value)
  oGRgenes.p = append(oGRgenes.p, gr[pos])
  names(oGRgenes.p)[i] = names(oGRLpeaks)[i]
}

# write gff file and make some plots
hist(oGRgenes.p$neg.log10.fisher.p.value, main='p value distribution', xlab='neg log10 fisher p.value')
temp = f_step6(f_step5(lCas, oGRgenes.p), main='Overlaps')

oGRgenes.p = oGRgenes.p[order(oGRgenes.p$neg.log10.fisher.p.value, decreasing = T)]

## create a gff and csv file for these 
r = oGRgenes.p
n = make.names(paste('genes_ranking_over_conserved_peaks', date(), 'gff'))
dir.create('Results_scripts/Gff_files', showWarnings = F)
export(r, con=paste('Results_scripts/Gff_files/', n, sep=''), format = 'gff3')
n = make.names(paste('genes_ranking_over_conserved_peaks', date(), 'csv'))
write.csv(as.data.frame(r), file=paste('Results_scripts/', n, sep=''))  


