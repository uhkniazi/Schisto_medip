# File: 05_RepEnrich_RNASeq_Analysis.R
# Auth: u.niazi@imperial.ac.uk
# DESC: Uses the summary of counts overlapping each repeat class in each sample by RepEnrich for the Treated and control libraries
#       and calculates proprtions by modelling the data as a multinomial distribution
# Date: 2/07/2015

source('00_header.R')

# import the data
dfFile = read.csv(file.choose(), header=T, stringsAsFactors=F)

mDat.import = matrix(NA, nrow = nrow(dfFile), ncol=9)
rownames(mDat.import) = dfFile$Sample

for (i in 1:nrow(dfFile)){
  f = read.csv(dfFile$File[i], sep='\t', header=F)
  # set column names
  if (i == 1) colnames(mDat.import) = f$V1
  mDat.import[i,] = f$V2
}

mDat = round(mDat.import)

###################################
## data normalization for library size
## equal distribution across each sample for a repeat class
mDat = t(log(mDat))
ivRM = rowMeans(mDat)
mDat = mDat - ivRM 
ivSF = apply(mDat, 2, median)
ivSF = exp(ivSF)

mDat = round(mDat.import)
#mDat = mDat/ivSF
mDat = sweep(mDat, MARGIN = 1, STATS = ivSF, FUN = '/')

#################################
## distribution across each repeat class for a sample
mDat.p = mDat/rowSums(mDat)
# take mean for each sample replicate
C = colMeans(mDat.p[1:3,])
S = colMeans(mDat.p[4:6,])
#######################
## function to simulate
f_sim_multi = function(prob){
  m = t(rmultinom(n = 1000, size = 1000, prob = prob))
  # convert to probability scale for plotting
  return(m/rowSums(m))
}

f_sim_ci = function(prob){
  m = f_sim_multi(prob)
  return(apply(m, 2, function(x) quantile(x, c(0.025, 0.975))))
}

#################################
## plots with error bars
mBar = rbind(C,S)
c = rainbow(nrow(mBar))
l2 = barplot(mBar, beside=T, las=2, col=c, main='Proportion of Reads Distributed Over Repeat Classes', ylim=c(0,1))
legend('topright', legend = rownames(mBar), fill=c)

# draw error bars
for (i in 1:nrow(mBar)){
  l = l2[i,]
  ## make error bars
  m = f_sim_ci(mBar[i,])
  segments(l, y0 = m[1,], l, y1 = m[2,], lwd=2)
  segments(l-0.1, y0 = m[1,], l+0.1, y1 = m[1,], lwd=2)
  segments(l-0.1, y0 = m[2,], l+0.1, y1 = m[2,], lwd=2)
}

# print confidence intervals
for (i in 1:nrow(mBar)){
  print(rownames(mBar)[i])
  print(round(f_sim_ci(mBar[i,]), 3))
}

## plot separately
mBar = mDat.p
c = rainbow(nrow(mBar))
l2 = barplot(mBar, beside=T, las=2, col=c, main='Proportion of Reads Distributed Over Repeat Classes', ylim=c(0,1))
legend('topright', legend = rownames(mBar), fill=c)

# draw error bars
for (i in 1:nrow(mBar)){
  l = l2[i,]
  ## make error bars
  m = f_sim_ci(mBar[i,])
  segments(l, y0 = m[1,], l, y1 = m[2,], lwd=2)
  segments(l-0.1, y0 = m[1,], l+0.1, y1 = m[1,], lwd=2)
  segments(l-0.1, y0 = m[2,], l+0.1, y1 = m[2,], lwd=2)
}

# print confidence intervals
for (i in 1:nrow(mBar)){
  print(rownames(mBar)[i])
  print(round(f_sim_ci(mBar[i,]), 3))
}

########### try with DESeq2

library(DESeq2)
mDat.d = t(round(mDat.import))
f = gl(2, 3, labels = c('C', 'S'))
dfDesign = data.frame(condition=f, row.names=colnames(mDat.d))

oDseq = DESeqDataSetFromMatrix(mDat.d, dfDesign, design = ~ condition)
oDseq = DESeq(oDseq)

res = results(oDseq)
