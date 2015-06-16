# File: 04_RepEnrich_MeDIP_Analysis.R
# Auth: u.niazi@imperial.ac.uk
# DESC: Uses the summary of counts overlapping each repeat class in each sample by RepEnrich for the MeDIP and control libraries
#       and calculates proprtions by modelling the data as a multinomial distribution
# Date: 12/05/2015


######################## 
## data import
dfDat = read.csv(file.choose(), header=T, row.names=1)

## data formatting
mDat = round(t(dfDat),0)
#mDat = mDat/rowSums(mDat)

###################################
## data normalization for library size
## equal distribution across each sample for a repeat class
mDat = t(log(mDat))
ivRM = rowMeans(mDat)
mDat = mDat - ivRM 
ivSF = apply(mDat, 2, median)
ivSF = exp(ivSF)

mDat = t(dfDat)
mDat = mDat/ivSF

#################################
## distribution across each repeat class for a sample
mDat.p = mDat/rowSums(mDat)
# take mean for each sample replicate
fme = colMeans(mDat.p[1:3,])
mme = colMeans(mDat.p[4:6,])
sme = colMeans(mDat.p[7:9,])

fun = colMeans(mDat.p[10:12,])
mun = colMeans(mDat.p[13:15,])
sun = colMeans(mDat.p[16:18,])

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
mBar = rbind(fme, mme, sme, fun, mun, sun)
c = rainbow(6)
l2 = barplot(mBar, beside=T, las=2, col=c, main='Proportion of Reads Distributed Over Repeat Classes', ylim=c(0,0.5))
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
