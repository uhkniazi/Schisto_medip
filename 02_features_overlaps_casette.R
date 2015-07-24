# File: 02_features_overlaps_casette.R
# Auth: u.niazi@imperial.ac.uk
# DESC: uses the casette object and medip peaks object created earlier, creates a gene casette and finds overlaps of
#       casettes and medip peaks. Reports statistics and plots
# Date: 01/05/2015


# source header file
source('00_header.R')

### data loading
p = paste('Choose R object file with pooled Medip peaks, look in Objects/ directoy')
print(p)
oGRpooled = f_LoadObject(file.choose())

p = paste('Choose features casette created earlier with create_features_from_gff.R')
print(p)
lFeatures = f_LoadObject(file.choose())

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
