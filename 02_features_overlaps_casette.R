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


