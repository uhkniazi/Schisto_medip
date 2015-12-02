# File: 02_features_overlaps_casette_2.R
# Auth: u.niazi@imperial.ac.uk
# DESC: uses the casette object and pooled medip peaks for each lifecycle stage created earlier, 
#       creates a gene casette and finds overlaps of casettes and medip peaks. Reports statistics and plots
# Date: 01/05/2015


# source header file
source('00_header.R')

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


# somules
f = f_step1(lCas, oGRLpooled$s)
s = f_step3(lCas, oGRLpooled$s[f])
s = f_step4(lCas, s)
s = f_step5(lCas, s)
# plot the results
mSom = f_step6(s, main='Somule')
somule = s$peaks


# female
f = f_step1(lCas, oGRLpooled$f)
s = f_step3(lCas, oGRLpooled$f[f])
s = f_step4(lCas, s)
s = f_step5(lCas, s)
# plot the results
mFem = f_step6(s, main='Female')
female = s$peaks

# male
f = f_step1(lCas, oGRLpooled$m)
s = f_step3(lCas, oGRLpooled$m[f])
s = f_step4(lCas, s)
s = f_step5(lCas, s)
# plot the results
mMale = f_step6(s, main='Male')
male = s$peaks

# plot them together
# get the median to plot
s = apply(mSom, 2, median)
f = apply(mFem, 2, median)
m = apply(mMale, 2, median)

mBar = rbind(s, f, m)
col = grey.colors(nrow(mBar))
l = barplot(mBar, beside=T, ylim=c(0,1), col=col, main='Probability of Peak over a Feature')
## draw error bars
f_barplot_errorbars = function(x.loc, y.loc, ...){
  segments(x.loc, y.loc[1], x.loc, y.loc[2], ...)
  segments(x.loc-0.1, y.loc[1], x.loc+0.1, y.loc[1], ...)
  segments(x.loc-0.1, y.loc[2], x.loc+0.1, y.loc[2], ...)
}
sapply(seq_along(1:ncol(mBar)), function(x) f_barplot_errorbars(l[1,x], quantile(mSom[,x], c(0.025, 0.975))))
sapply(seq_along(1:ncol(mBar)), function(x) f_barplot_errorbars(l[2,x], quantile(mFem[,x], c(0.025, 0.975))))
sapply(seq_along(1:ncol(mBar)), function(x) f_barplot_errorbars(l[3,x], quantile(mMale[,x], c(0.025, 0.975))))

legend('topright', legend = c(rownames(mBar)), fill=col)


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










## data processing
# remove the medip peaks from class 'none'
f = which(oGRpooled$groups.lab == 'none')
oGRpooled = oGRpooled[-f]
# create a matrix to hold the matches
mFeatures = matrix(NA, nrow = length(lFeatures$fst.exon), ncol = 5)
colnames(mFeatures) = c('upstream', 'fst.exon', 'exons.others', 'introns', 'downstream')
rownames(mFeatures) = names(lFeatures$fst.exon)
cn = colnames(mFeatures)
rn = rownames(mFeatures)
# fill the data
for(i in 1:ncol(mFeatures)){
  # fill data one column at a time
  mFeatures[,i] = countOverlaps(lFeatures[[cn[i]]][rn], oGRpooled)
}

# summarize the data
mBar = colSums(mFeatures)
print(mBar)
mBar = mBar / sum(mBar)
barplot(mBar)

f = rowSums(mFeatures)
f = which(f == 0)

mFeatures.sub = mFeatures[-f,]

write.csv(data.frame(mFeatures.sub), file='Reports/casette_genes.csv')

### create groups of this data based on pattern of methylation
# easy way to cut is using a hclust function
mOverlap.2.s = t(apply(mFeatures.sub, 1, as.logical))
dm = dist(mOverlap.2.s, method = 'binary')
hc = hclust(dm)
#plot(hc, labels=F, hang=0)
cp = cutree(hc, h = 0.18)
# get a new matrix based on these cutpoint members
mGroups = NULL # this is a matrix of all possible combinations/classes
for (i in 1:length(unique(cp))){
  i2 = which(cp == i)
  v = NULL
  if(length(i2) == 1) v = mOverlap.2.s[i2,] else v= as.logical(colSums(mOverlap.2.s[i2,]))
  mGroups = rbind(mGroups, v)
}
colnames(mGroups) = colnames(mFeatures.sub)
n = c('i', 'ds', 'us', 'i.ds', 'ex.i.ds', 'ex.i', 'us.i', 
      'us.fex.i', 'us.ds', 'us.fex.i.ds', 'us.fex', 'ex.ds',
      'ex', 'us.ex.i', 'fex.i', 'fex', 'us.fex.ex.i', 'us.ex.ds',
      'us.fex.ex.i.ds', 'fex.ex.i', 'us.fex.ds', 'us.ex', 'fex.ex.i.ds',
      'us.ex.i.ds', 'fex.ds')
rownames(mGroups) = n
print('Possible groups in the data')
print(mGroups)
# assign a group to each member of the full matrix
groups = rep(NA, length.out=nrow(mOverlap.2.s))
for (i in 1:length(unique(cp))){
  f = apply(mOverlap.2.s, 1, function(x) all(x == mGroups[i,]))
  groups[f] = rownames(mGroups)[i]
}
#groups.lab = factor(groups, labels = c('sfm', 'sf', 'sm', 's', 'f', 'm', 'fm'))
groups.lab = factor(groups)
print('Group count'); print(table(groups))
print('Group count with labels'); print(table(groups.lab))
# final data frame
dfCasetteGroups = data.frame(groups.lab, mFeatures.sub)

## which ones overlap with particular lifecycle stages
mFeatures = matrix(NA, nrow = length(lFeatures$fst.exon), ncol = 5)
colnames(mFeatures) = c('upstream', 'fst.exon', 'exons.others', 'introns', 'downstream')
rownames(mFeatures) = names(lFeatures$fst.exon)
cn = colnames(mFeatures)
rn = rownames(mFeatures)
i2 = grep('s', oGRpooled$groups.lab)
# fill the data
for(i in 1:ncol(mFeatures)){
  # fill data one column at a time
  mFeatures[,i] = countOverlaps(lFeatures[[cn[i]]][rn], oGRpooled[i2])
}
f = rowSums(mFeatures)
f = which(f == 0)
mFeatures.sub.somule = mFeatures[-f,]
i = rownames(dfCasetteGroups) %in% rownames(mFeatures.sub.somule)
dfCasetteGroups$somule = i

# repeat for female
mFeatures = matrix(NA, nrow = length(lFeatures$fst.exon), ncol = 5)
colnames(mFeatures) = c('upstream', 'fst.exon', 'exons.others', 'introns', 'downstream')
rownames(mFeatures) = names(lFeatures$fst.exon)
cn = colnames(mFeatures)
rn = rownames(mFeatures)
i2 = grep('f', oGRpooled$groups.lab)
# fill the data
for(i in 1:ncol(mFeatures)){
  # fill data one column at a time
  mFeatures[,i] = countOverlaps(lFeatures[[cn[i]]][rn], oGRpooled[i2])
}
f = rowSums(mFeatures)
f = which(f == 0)
mFeatures.sub.female = mFeatures[-f,]
i = rownames(dfCasetteGroups) %in% rownames(mFeatures.sub.female)
dfCasetteGroups$female = i

# repeat for male
mFeatures = matrix(NA, nrow = length(lFeatures$fst.exon), ncol = 5)
colnames(mFeatures) = c('upstream', 'fst.exon', 'exons.others', 'introns', 'downstream')
rownames(mFeatures) = names(lFeatures$fst.exon)
cn = colnames(mFeatures)
rn = rownames(mFeatures)
i2 = grep('m', oGRpooled$groups.lab)
# fill the data
for(i in 1:ncol(mFeatures)){
  # fill data one column at a time
  mFeatures[,i] = countOverlaps(lFeatures[[cn[i]]][rn], oGRpooled[i2])
}
f = rowSums(mFeatures)
f = which(f == 0)
mFeatures.sub.male = mFeatures[-f,]
i = rownames(dfCasetteGroups) %in% rownames(mFeatures.sub.male)
dfCasetteGroups$male = i

write.csv(dfCasetteGroups, file='Reports/casette_genes_with_groups.csv')
save(dfCasetteGroups, file='Objects/dfCasetteGroups.rds')
