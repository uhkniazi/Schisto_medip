# File: 04_RepEnrich_MeDIP_Analysis.R
# Auth: u.niazi@imperial.ac.uk
# DESC: Uses the summary of counts overlapping each repeat class in each sample by RepEnrich for the MeDIP and control libraries
#       and calculates proprtions by modelling the data as a multinomial distribution
# Date: 12/05/2015

source('00_header.R')

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
#mDat = mDat/ivSF
mDat = sweep(mDat, MARGIN = 1, STATS = ivSF, FUN = '/')


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
l2 = barplot(mBar, beside=T, las=2, col=c, main='Proportion of Reads Distributed Over Repeat Classes', ylim=c(0,0.6))
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


#########################################################################
#### repeat analysis with real medip data

## load the medip data
oGRpooled = f_LoadObject(file.choose())
oGRpooled = oGRpooled[oGRpooled$groups.lab != 'none']
oGRpooled = oGRpooled[seqnames(oGRpooled) %in% gcvChromosomes]
# load the repeats object 
oGRrep = f_LoadObject(file.choose())
strand(oGRrep) = '*'
# rename all the W chromosome to ZW
r1 = oGRrep[seqnames(oGRrep) %in% gcvChromosomes]
r2 = oGRrep[seqnames(oGRrep) == 'Schisto_mansoni.Chr_W']
r1 = GRanges(as.character(seqnames(r1)), IRanges(start(r1), end(r1)), strand(r1), DataFrame(mcols(r1)))
r2 = GRanges('Schisto_mansoni.Chr_ZW', IRanges(start(r2), end(r2)), strand(r2), DataFrame(mcols(r2)))
oGRrep = c(r1, r2)

f = oGRrep$mcols.class
oGRLrep = split(oGRrep, f)

# repeat sequences reads ratio from RepEnrich
ivRep = colMeans(mDat.p[10:18,])
# count overlaps
mMed = sapply(oGRLrep, function(x) factor(overlapsAny(x, oGRpooled),levels = c('TRUE', 'FALSE')))
mMed = sapply(mMed, table)

ivMed = mMed['TRUE',] / rowSums(mMed)['TRUE']
ivMed = ivMed[names(ivRep)]
# mMed = mMed[,names(ivRep)]
#### plot the data
mBar = rbind(ivMed, ivRep)
rownames(mBar) = c('MeDIP', 'Unmethylated')
c = rainbow(2)

l2 = barplot(mBar, beside=T, las=2, main='Proportion of Reads Distributed Over Repeat Classes', col=c, ylim=c(0, 0.6),
             ylab='fraction')
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

# compare the repeat classes to 
# the background i.e. unmethylated classes
# total trials
n = rowSums(mMed)['TRUE']
# under the null hypothesis the expected proportion should look like
df = data.frame(Expected.prop=round(ivRep, 3), Expected.count=round(n*ivRep, 0))
df = data.frame(Observed.prop=round(ivMed, 3), Observed.count=mMed['TRUE', names(ivRep)])

p.vals = rep(NA, length=(length(ivRep)))
names(p.vals) = names(ivRep)
for (i in 1:length(ivRep))
{
  nm = names(ivRep)[i]
  x = mMed['TRUE', nm]
  p.vals[i] = binom.test(x, n, p = ivRep[nm])$p.value
  print(paste(names(ivRep)[i], signif(p.vals[i], 3)))
}

p.adj = p.adjust(p.vals, method = 'BH')
as.data.frame(print(signif(p.adj, 3)))
# # fisher odds ratio test
# p.vals = rep(NA, length=(length(ivRep)))
# names(p.vals) = names(ivRep)
# for (i in 1:length(ivRep)){
#   o1 = c(mMed['TRUE', i], n-mMed['TRUE',i])
#   o2 = round(ivRep[i] * n)
#   o2 = c(o2, n-o2)
#   m = matrix(c(o1, o2), nrow=2)
#   p.vals[i] = fisher.test(m, simulate.p.value = T)$p.value
# }


#######################################################################
##### breakdown into lifecycle stages

## female
# repeat sequences reads ratio from RepEnrich
ivRep = fun
# count overlaps
f = which(oGRpooled$groups.lab %in% c('f', 'sfm', 'fm', 'sf'))
mMed = sapply(oGRLrep, function(x) factor(overlapsAny(x, oGRpooled[f]),levels = c('TRUE', 'FALSE')))
mMed = sapply(mMed, table)

ivMed = mMed['TRUE',] / rowSums(mMed)['TRUE']
ivMed = ivMed[names(ivRep)]

#### plot the data
mBar = rbind(ivMed, ivRep)
rownames(mBar) = c('MeDIP.F', 'fun')
c = rainbow(2)

l2 = barplot(mBar, beside=T, las=2, main='Proportion of Reads Distributed Over Repeat Classes', col=c, ylim=c(0, 0.6),
             ylab='fraction')
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


## male
# repeat sequences reads ratio from RepEnrich
ivRep = mun
# count overlaps
f = which(oGRpooled$groups.lab %in% c('m', 'sfm', 'fm', 'sm'))
mMed = sapply(oGRLrep, function(x) factor(overlapsAny(x, oGRpooled[f]),levels = c('TRUE', 'FALSE')))
mMed = sapply(mMed, table)

ivMed = mMed['TRUE',] / rowSums(mMed)['TRUE']
ivMed = ivMed[names(ivRep)]

#### plot the data
mBar = rbind(ivMed, ivRep)
rownames(mBar) = c('MeDIP.M', 'mun')
c = rainbow(2)

l2 = barplot(mBar, beside=T, las=2, main='Proportion of Reads Distributed Over Repeat Classes', col=c, ylim=c(0, 0.6),
             ylab='fraction')
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

## somule
# repeat sequences reads ratio from RepEnrich
ivRep = sun
# count overlaps
f = which(oGRpooled$groups.lab %in% c('s', 'sfm', 'sm', 'sf'))
mMed = sapply(oGRLrep, function(x) factor(overlapsAny(x, oGRpooled[f]),levels = c('TRUE', 'FALSE')))
mMed = sapply(mMed, table)

ivMed = mMed['TRUE',] / rowSums(mMed)['TRUE']
ivMed = ivMed[names(ivRep)]

#### plot the data
mBar = rbind(ivMed, ivRep)
rownames(mBar) = c('MeDIP.S', 'sun')
c = rainbow(2)

l2 = barplot(mBar, beside=T, las=2, main='Proportion of Reads Distributed Over Repeat Classes', col=c, ylim=c(0, 0.6),
             ylab='fraction')
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

######################################################################
### life cycle stages together with unmethylated
# repeat sequences reads ratio from RepEnrich
ivRep = colMeans(mDat.p[10:18,])

f = which(oGRpooled$groups.lab %in% c('f', 'sfm', 'fm', 'sf'))
mMed = sapply(oGRLrep, function(x) factor(overlapsAny(x, oGRpooled[f]),levels = c('TRUE', 'FALSE')))
mMed = sapply(mMed, table)
ivMed = mMed['TRUE',] / rowSums(mMed)['TRUE']
ivMed = ivMed[names(ivRep)]
ivMed.f = ivMed

f = which(oGRpooled$groups.lab %in% c('m', 'sfm', 'fm', 'sm'))
mMed = sapply(oGRLrep, function(x) factor(overlapsAny(x, oGRpooled[f]),levels = c('TRUE', 'FALSE')))
mMed = sapply(mMed, table)
ivMed = mMed['TRUE',] / rowSums(mMed)['TRUE']
ivMed = ivMed[names(ivRep)]
ivMed.m = ivMed

f = which(oGRpooled$groups.lab %in% c('s', 'sfm', 'sm', 'sf'))
mMed = sapply(oGRLrep, function(x) factor(overlapsAny(x, oGRpooled[f]),levels = c('TRUE', 'FALSE')))
mMed = sapply(mMed, table)
ivMed = mMed['TRUE',] / rowSums(mMed)['TRUE']
ivMed = ivMed[names(ivRep)]
ivMed.s = ivMed

#### plot the data
mBar = rbind(ivMed.s, ivMed.f, ivMed.m, ivRep)
rownames(mBar) = c('MeDIP.S', 'MeDIP.F', 'MeDIP.M', 'Unmethylated')
c = rainbow(nrow(mBar))

l2 = barplot(mBar, beside=T, las=2, main='Proportion of Reads Distributed Over Repeat Classes', col=c, ylim=c(0, 0.6),
             ylab='fraction')
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

