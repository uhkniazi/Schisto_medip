# File: 03_medip_overlap_genome.R
# Auth: u.niazi@imperial.ac.uk
# DESC: uses the medip data and features (including repeats) objects created earlier and finds overlaps
# Date: 09/05/2015


# source header file
source('00_header.R')

### data loading
p = paste('Choose R object file with pooled Medip peaks, look in Objects/ directoy')
print(p)
oGRpooled = f_LoadObject(file.choose())

p = paste('Choose features object created earlier with create_features_from_gff.R')
print(p)
lFeatures = f_LoadObject(file.choose())

#########################################################
## data processing
# remove the medip peaks from class 'none'
f = which(oGRpooled$groups.lab == 'none')
oGRpooled = oGRpooled[-f]

## count over features granges 
lFeatures$gene = NULL
lFeatures$desc = NULL
lOverlaps = sapply(lFeatures, function(x) overlapsAny(x, oGRpooled))
mDat = sapply(lOverlaps, table)

rownames(mDat) = c('False', 'True')

# calculate confidence intervals/errors by simulating the data
ivProb.medip.all = mDat['True',] / rowSums(mDat)['True']

# simulate data as a sample from a multinomial distribution
# P(y | Theta)
mDat.medip.all = mDat['True',]
mSim.medip.all = t(rmultinom(n = 1000, size = 1000, prob = ivProb.medip.all))
# convert to probability scale for plotting
mSim.medip.all = mSim.medip.all / rowSums(mSim.medip.all)
# look at data distribution
l = barplot(colMeans(mSim.medip.all), ylim = c(0, 0.8), main='Distribution of all peaks')
# calculate quantiles for drawing
m = apply(mSim.medip.all, 2, function(x) quantile(x, c(0.025, 0.975)))
segments(l, y0 = m[1,], l, y1 = m[2,], lwd=2)
segments(l-0.1, y0 = m[1,], l+0.1, y1 = m[1,], lwd=2)
segments(l-0.1, y0 = m[2,], l+0.1, y1 = m[2,], lwd=2)

#######################################################
### repeat this analysis for all classes of medip peaks
lCI = vector('list', length = 7)
lProb = vector('list', length = 7)
n = c("sf", "s", "m", "sfm", "f", "sm", "fm")
names(lCI) = n
names(lProb) = n
cvNames = names(lCI)
# loop through each class of medip peak and calculate the proportions
for (i in 1:length(lCI)){
  # select pooled peaks from particular class
  lOverlaps = sapply(lFeatures, function(x) overlapsAny(x, oGRpooled[oGRpooled$groups.lab == cvNames[i]]))
  mDat = sapply(lOverlaps, table)
  rownames(mDat) = c('False', 'True')
  # calculate confidence intervals/errors by simulating the data
  ivProb.medip = mDat['True',] / rowSums(mDat)['True']
  lProb[[cvNames[i]]] = ivProb.medip  
  # simulate data as a sample from a multinomial distribution
  # P(y | Theta)
  mDat.medip = mDat['True',]
  mSim.medip = t(rmultinom(n = 1000, size = 1000, prob = ivProb.medip))
  # convert to probability scale for plotting and calculate 95% CI
  mSim.medip = mSim.medip / rowSums(mSim.medip)
  # calculate quantiles for drawing
  lCI[[cvNames[i]]] = apply(mSim.medip, 2, function(x) quantile(x, c(0.025, 0.975)))
}

## collapse the list into a matrix
col = rainbow(length(cvNames))
mBar = sapply(lProb, function(x) x)
mBar = t(mBar)

l2 = barplot(mBar, beside=T, ylim=c(0, 0.8), col=col)
legend('topleft', legend = rownames(mBar), fill = col)

# draw error bars
for(i in 1:nrow(l2)){
  l = l2[i,]
  m = lCI[[cvNames[i]]]
  segments(l, y0 = m[1,], l, y1 = m[2,], lwd=2)
  segments(l-0.1, y0 = m[1,], l+0.1, y1 = m[1,], lwd=2)
  segments(l-0.1, y0 = m[2,], l+0.1, y1 = m[2,], lwd=2)
}

sapply(seq_along(lCI), function(x) {
  print(names(lCI)[x])
  print(round(lCI[[x]], 2))
})
























