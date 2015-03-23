# File: 01_dismiss_to_gff.R
# Auth: umn@aber.ac.uk
# DESC: requires a csv file with 2 columns - one representing sample names and second representing
#       file paths to the .rds objects created by dismiss. it divides the data into groups depending on 
#       presence/absence in each sample i.e. Somule, Female and Male with possible values of True and false
#       for each condition - thus giving 2^3 - 1 possible groups - which would be 7 in our case.
#       the output is a gff file with location, strand and group of each peak.
# Date: 9/12/2014

### input variables setting
# filter size to use later on, set to 0 or 1
# 0 means a peak should be present in at least 1 replicate
# 1 means a peak should be present in at least 2 replicates
iFilterSize = 2

# get the csv file with sample names and file paths 
# this will be data output from dismiss
# example of first 2 lines of the csv file, the first line is the header
# sample  path
# s1  /mnt/windrv/Shared/Data/Project_Medip/Data_peak_calling/Macs2.s1/oGRdismiss_Fri.Nov.28.10.23.02.2014.rds
print('Select CSV file with sample names and dismiss output .rds files')
setwd(pPeaks)
# load the csv file with samples
dat.file = file.choose()
dfSamples = read.csv(dat.file, header=T)

### processing steps
# save old plotting parameters
p.old = par()
lSamples = vector('list', length=nrow(dfSamples))
names(lSamples) = dfSamples$sample
n = names(lSamples)
for (i in 1:nrow(dfSamples)) lSamples[[n[i]]] = f_LoadObject(as.character(dfSamples[i,'path']))
setwd(pWD)

# make granges list of medip-seq peaks data
oGRLsamples = GRangesList(lSamples)
names(oGRLsamples) = names(lSamples)

# how many peaks in each sample
print('The number of peaks in each sample')
print((table(strand(oGRLsamples))))
# plot number of peaks in each sample
# what is the distribution of plus minus and double stranded peaks
t = (table(strand(oGRLsamples)))
c = grey.colors(3)
yl = ceiling((max(t)/1000))*1000
barplot(t(t), beside=T, ylim=c(0, yl), col = c,
        xlab='samples', ylab='count', main='number of peaks in each sample')
legend('topleft', legend = colnames(t), fill = c)

# the average sizes of peaks should be around 300-500 as that is the resolution of medip-seq
# sizes of peaks
w = width(oGRLsamples)
print(paste('average peak sizes'))
print(sapply(w, summary))
par(mfrow=c(3,3))
n = names(w)
for (i in 1:length(w)) {
  t = paste('peak width distribution in', n[i])
  w2 = w[[i]]
  w2 = w2[w2 < quantile(w2, 0.90)]
  hist(w2, main=t, xlab='width', ylab='freq')
}
par(p.old)

# peak distribution and p.values
n = names(oGRLsamples)
par(mfrow=c(1,2))
for (i in 1:length(oGRLsamples)){
  p = oGRLsamples[[i]]$neg.log10.qval
  s = as.factor(strand(oGRLsamples[[i]]))
  boxplot(p ~ s, main=paste('qvalues distribution by strands in', n[i]), 
          xlab='strand', ylab='neg log 10 qvalue')
  stripchart(p ~ s, main=paste('qvalues distribution by strands in', n[i]), 
             ylab='strand', xlab='neg log 10 qvalue', method='jitter')
}

# how many strands overlap
df = f_dfGetMatchingStrandsFromGRanges(oGRLsamples[[1]], oGRLsamples[[2]])
t = table(s1=df$s.query, s2=df$s.subject)
print(t)
print(paste(sum(rowSums(t)), 'out of', length(oGRLsamples[[1]]) ))

df = f_dfGetMatchingStrandsFromGRanges(oGRLsamples[[3]], oGRLsamples[[2]])
t = table(s3=df$s.query, s2=df$s.subject)
print(t)
print(paste(sum(rowSums(t)), 'out of', length(oGRLsamples[[3]]) ))

df = f_dfGetMatchingStrandsFromGRanges(oGRLsamples[['f3']], oGRLsamples[['f2']])
t = table(f3=df$s.query, f2=df$s.subject)
print(t)
print(paste(sum(rowSums(t)), 'out of', length(oGRLsamples[['f3']]) ))

df = f_dfGetMatchingStrandsFromGRanges(oGRLsamples[['f4']], oGRLsamples[['f2']])
t = table(f4=df$s.query, f2=df$s.subject)
print(t)
print(paste(sum(rowSums(t)), 'out of', length(oGRLsamples[['f4']]) ))

df = f_dfGetMatchingStrandsFromGRanges(oGRLsamples[['m1']], oGRLsamples[['m3']])
t = table(m1=df$s.query, m3=df$s.subject)
print(t)
print(paste(sum(rowSums(t)), 'out of', length(oGRLsamples[['m1']]) ))

df = f_dfGetMatchingStrandsFromGRanges(oGRLsamples[['m2']], oGRLsamples[['m3']])
t = table(m2=df$s.query, m3=df$s.subject)
print(t)
print(paste(sum(rowSums(t)), 'out of', length(oGRLsamples[['m2']]) ))

## classify data into groups
## how many unique classes
gr = unlist(oGRLsamples)
print(paste('Total number of pooled peaks',length(gr)))
gr = reduce(gr)
print(paste('Total number of unique peaks',length(gr)))
# create a binary matrix with columns as samples and rows as methylated regions
mOverlap = matrix(NA, nrow = length(gr), ncol = length(oGRLsamples))
rownames(mOverlap) = 1:length(gr)
colnames(mOverlap) = names(oGRLsamples)
# check which vector (column vector) overlaps in which samples
# fill each column with overlap vector
n = names(oGRLsamples)
for (i in 1:length(n)){
  v = overlapsAny(gr, oGRLsamples[[n[i]]])
  mOverlap[,n[i]] = v
}

# a peak should be TRUE in 1 or 2 replicates at least - set iFilterSize variable
# adding the row values for replicates e.g. s1, s2 and s3 will give a count of 0 to 3
s = mOverlap[,'s1'] + mOverlap[,'s2'] + mOverlap[,'s3']
f = mOverlap[,'f2'] + mOverlap[,'f3'] + mOverlap[,'f4']
m = mOverlap[,'m1'] + mOverlap[,'m2'] + mOverlap[,'m3']
# convert to a boolean vector depending on cutoff point
s = s > iFilterSize
f = f > iFilterSize
m = m > iFilterSize

# join the 3 vectors into a matrix N X 3 
mOverlap.2 = cbind(s, f, m)
rownames(mOverlap.2) = 1:length(gr)

# not using this filter anymore
# # merge the sample columns i.e. if e.g. for somules it is true in any of the 3 
# # samples then use it as true or else false by using an OR function
# s = mOverlap[,'s1'] | mOverlap[,'s2'] | mOverlap[,'s3']
# f = mOverlap[,'f2'] | mOverlap[,'f3'] | mOverlap[,'f4']
# m = mOverlap[,'m1'] | mOverlap[,'m2'] | mOverlap[,'m3']
# 
# mOverlap.2 = cbind(s, f, m)
# rownames(mOverlap.2) = 1:length(gr)
par(p.old)
# take a sample of the rows to create classes using dist and hclust
# easy way to cut is using a hclust function
set.seed(123)
x = sample(1:length(gr), size = 1000, replace = F)
mOverlap.2.s = mOverlap.2[x,]

dm = dist(mOverlap.2.s, method = 'binary')
hc = hclust(dm)
#plot(hc, labels=F, hang=0)
cp = cutree(hc, h = 0.2) # or cut it as 7 as there are 2^3 -1 possible classes
# get a new matrix based on these cutpoint members
mGroups = NULL # this is a matrix of all possible combinations/classes
for (i in 1:length(unique(cp))){
  v = colSums(mOverlap.2.s[cp == i, ])
  mGroups = rbind(mGroups, as.logical(v))
}
colnames(mGroups) = colnames(mOverlap.2)
print('Possible groups in the data')
print(mGroups)
# assign a group to each member of the full matrix
groups = rep(NA, length.out=nrow(mOverlap.2))
for (i in 1:length(unique(cp))){
  # vectorized match by using a single &
  f = (mOverlap.2[,1] == mGroups[i,1]) & (mOverlap.2[,2] == mGroups[i,2]) & (mOverlap.2[,3] == mGroups[i,3])
  groups[f] = i
}
#groups.lab = factor(groups, labels = c('sfm', 'sf', 'sm', 's', 'f', 'm', 'fm'))
groups.lab = factor(groups, labels = c('none', 'sfm', 'sm', 'm', 'fm', 'f', 'sf', 's'))
print('Group count'); print(table(groups))
print('Group count with labels'); print(table(groups.lab))
# final data frame
dfPeakGroups = data.frame(groups.lab, groups, mOverlap.2)

oGRpooled = gr
mcols(oGRpooled) = dfPeakGroups
md = list(groups=mGroups, desc='groups created by 01_dismiss_to_gff.R using iFilter = 2', date=paste(date()))
metadata(oGRpooled) = md
n = make.names(paste('oGRpooled.medip.peaks', date(), '.rds'))
dir.create('Objects', showWarnings = F)
save(oGRpooled, file=paste('Objects/', n, sep=''))

## create a gff of the pooled peaks and their groups
n = make.names(paste('oGRpooled.medip.peaks', date(), '.gff'))
dir.create('Results_scripts/Gff_files', showWarnings = F)
export(oGRpooled, con=paste('Results_scripts/Gff_files/', n, sep=''), format = 'gff3')
