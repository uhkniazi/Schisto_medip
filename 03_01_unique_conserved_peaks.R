# File: 03_01_unique_conserved_peaks.R
# Auth: u.niazi@imperial.ac.uk
# DESC: uses the pooled medip peaks from each lifecycle stage, created earlier using 01_02_combine_replicates.R and gets the unique
#       and conserved peaks.
# Date: 06/01/2016


# source header file
source('00_header.R')

### data loading
## pooled lifecycle stages with fisher p.values
oGRLpooled = f_LoadObject(file.choose())

# sanity check
sapply(oGRLpooled, length)
metadata(oGRLpooled)


### steps for analysis
# step 1 - get pooled peaks by union operation
# step 2 - matrix of overlaps
# step 3 - group the peaks
# step 4 - map the pooled peaks to lifecycle stage
# step 5 - match conserved peaks between the 3 samples and assign p-values


### functions used in script
f_step1 = function(p){
  u1 = union(p[[1]], p[[2]])
  u2 = union(u1, p[[3]])
  return(u2)
}

f_step2 = function(p, a){
  # match the strands
  m = matrix(FALSE, nrow=length(p), ncol = length(a))
  colnames(m) = names(a)
  # match the strands and positions
  for(i in 1:ncol(m)) {
    df = f_dfGetMatchingStrandsFromGRanges(p, a[[i]])
    # select only those where strands match
    f = df$s.query == df$s.subject
    pos = df$queryHits[f]
    # due to the width of the query peak, it may match more than one subject
    m[pos,i] = TRUE
  }
  return(m)
}

f_step3 = function(mat){
  s = c(T, F, F)
  f = c(F, T, F)
  m = c(F, F, T)
  cons = c(T, T, T)
  none = rep('none', times=nrow(mat))
  b = apply(mat, 1, function(x) all(x == s))
  none[b] = 's'
  b = apply(mat, 1, function(x) all(x == f))
  none[b] = 'f'
  b = apply(mat, 1, function(x) all(x == m))
  none[b] = 'm'
  b = apply(mat, 1, function(x) all(x == cons))
  none[b] = 'cons'
  return(none)
}

f_step4 = function(p, s){
  df = f_dfGetMatchingStrandsFromGRanges(s, p)
  # select only those where strands match
  f = df$s.query == df$s.subject
  df = df[f,]
  # remove duplicates i.e. peaks matching multiple positions in pooled data
  pos = df$queryHits[!duplicated(df$queryHits)]
  s = s[pos]
  return(s)
}

f_step5 = function(p, s){
  s = sapply(s, function(x) f_step4(p, x))
  # some wide pooled peaks overlap the individual samples, but the samples
  # themselves do not overlap, due to the union function used at step1
  # remove those peaks
  df.f = f_dfGetMatchingStrandsFromGRanges(s[['s']], s[['f']])
  df.m = f_dfGetMatchingStrandsFromGRanges(s[['s']], s[['m']])
  # remove duplicates and non matching strands
  # select only those where strands match
  f = df.f$s.query == df.f$s.subject
  df.f = df.f[f,]
  f = df.m$s.query == df.m$s.subject
  df.m = df.m[f,]
  # remove duplicates
  f = !duplicated(df.f$queryHits)
  df.f = df.f[f,]
  f = !duplicated(df.f$subjectHits)
  df.f = df.f[f,]
  
  f = !duplicated(df.m$queryHits)
  df.m = df.m[f,]
  f = !duplicated(df.m$subjectHits)
  df.m = df.m[f,]
  
  # match the 2 dataframes
  i = match(df.f$queryHits, df.m$queryHits)
  df = data.frame(cbind(s1=df.f$queryHits, s2=df.m$queryHits[i], f=df.f$subjectHits, m=df.m$subjectHits[i]))
  df = na.omit(df)
  som = s[['s']][df$s1]
  fem = s[['f']][df$f]
  mal = s[['m']][df$m]
  # calculate unified p-values
  p.val = rep(NA, times=nrow(df))
  for (i in seq_along(p.val)){
    p = c(som$q.value.1[i], som$q.value.2[i], fem$q.value.1[i], fem$q.value.2[i], mal$q.value.1[i], mal$q.value.2[i])
    p.val[i] = f_fishersMethod(p)
  }
  ret = union(som, fem)
  ret = union(ret, mal)
  ret$fisher.p.value = p.val
  ret$neg.log10.fisher.p.value = -1 * log10(p.val)
  return(ret)
}

oGRpooled = f_step1(oGRLpooled)
mPeaks = f_step2(oGRpooled, oGRLpooled)
fGroups = f_step3(mPeaks)
oGRpooled$fGroups = fGroups
oGRsomule = f_step4(oGRpooled[fGroups == 's'], oGRLpooled[['s']])
oGRfemale = f_step4(oGRpooled[fGroups == 'f'], oGRLpooled[['f']])
oGRmale = f_step4(oGRpooled[fGroups == 'm'], oGRLpooled[['m']])
oGRconserved = f_step5(oGRpooled[fGroups == 'cons'], oGRLpooled)

oGRLunique = GRangesList(s=oGRsomule, f=oGRfemale, m=oGRmale)


### report some basic stats
sapply(oGRLunique, length)
sapply(oGRLunique, function(x) table(strand(x)))
sapply(oGRLunique, function(x){
  df = as.data.frame(table(seqnames(x)))
  df = df[df$Freq > 10,]
  print(df)
})

p.old = par(mfrow=c(2,2))
sapply(seq_along(oGRLunique), function(x){
  p = oGRLunique[[x]]$neg.log10.fisher.p.value
  hist(p, main=names(oGRLunique)[x], xlab='Negative Log10 Fisher P-value')
})

length(oGRconserved)
table(strand(oGRconserved))
df = as.data.frame(table(seqnames(oGRconserved)))
print(df[df$Freq > 50,])
par(p.old)
hist(oGRconserved$neg.log10.fisher.p.value, main='Conserved Peaks', xlab='Negative Log10 Fisher P-value')


### save data
md = list(desc='Unique peaks for each lifecycle stage with fisher p.values', date=paste(date()))
metadata(oGRLunique) = md
n = make.names(paste('oGRLunique', date(), '.rds'))
dir.create('Objects', showWarnings = F)
save(oGRLunique, file=paste('Objects/', n, sep=''))

## create a gff for each sample
for (i in 1:length(oGRLunique)){
  r = oGRLunique[[i]]
  mc = as.data.frame(mcols(r))
  mc = data.frame(apply(mc, 2, function(x) signif(x, 3)))
  mcols(r) = mc
  n = make.names(paste('unique', names(oGRLunique)[i], date(), 'gff'))
  dir.create('Results_scripts/Gff_files', showWarnings = F)
  export(r, con=paste('Results_scripts/Gff_files/', n, sep=''), format = 'gff3')
  n = make.names(paste('unique', names(oGRLunique)[i], date(), 'csv'))
  write.csv(as.data.frame(r), file=paste('Results_scripts/', n, sep=''))  
}


## conserved peaks
md = list(desc='Conserved peaks common between lifecycle stage with fisher p.values', date=paste(date()))
metadata(oGRconserved) = md
n = make.names(paste('oGRconserved', date(), '.rds'))
dir.create('Objects', showWarnings = F)
save(oGRconserved, file=paste('Objects/', n, sep=''))

## create a gff for each sample
r = oGRconserved
mc = as.data.frame(mcols(r))
mc = data.frame(apply(mc, 2, function(x) signif(x, 3)))
mcols(r) = mc
n = make.names(paste('conserved', 'peaks', date(), 'gff'))
dir.create('Results_scripts/Gff_files', showWarnings = F)
export(r, con=paste('Results_scripts/Gff_files/', n, sep=''), format = 'gff3')
n = make.names(paste('conserved', 'peaks', date(), 'csv'))
write.csv(as.data.frame(r), file=paste('Results_scripts/', n, sep=''))





