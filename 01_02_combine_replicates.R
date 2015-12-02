# File: 01_02_combine_replicates.R
# Auth: uniazi@imperial.ac.uk
# DESC: Takes the oGRLsamples.fisher created in script 01_dismiss_to_gff.R and combines lifecycle replicates
#       into one representative for the lifecycle stage
# Date: 17/11/15



source('00_header.R')

# load the object created earlier
oGRLsamples.fisher = f_LoadObject(file.choose())

# sanity check if the correct object loaded
names(oGRLsamples.fisher)
metadata(oGRLsamples.fisher)
sapply(oGRLsamples.fisher, length)

## for combining the data certain parameters need to be addressed
## details can be seen here
## https://www.evernote.com/shard/s288/nl/38698211/19bc74ad-64c1-4977-bf68-5ac5b734dd5c
## step 1 - select peaks that are overlapping & have the same strands
## step 2 - intersect those peaks
## step 3 - carry forward the metadata fisher.p.value and neg.log10.fisher.p.value

# returns the boolean vector for regions that are true
f_step1 = function(g1, g2){
  df = f_dfGetMatchingStrandsFromGRanges(g1, g2)
  # very wide peaks can be duplicated, so remove those
  d = !duplicated(df$queryHits)
  df = df[d,]
  # get those regions where step 1 TRUE
  f = df$s.query == df$s.subject
  return(f)
}

# returns the intersection of the selected regions
f_step2 = function(g1, g2){
  return(intersect(g1, g2))
}

# get the metadata
f_step3 = function(g1, g2){
  mc = as.data.frame(mcols(g1))
  mc = mc[,c('fisher.p.value', 'neg.log10.fisher.p.value')]
  return(mc)
}

## somules
f.som = f_step1(oGRLsamples.fisher$s2, oGRLsamples.fisher$s3)
s.com = f_step2(oGRLsamples.fisher$s2[f.som], oGRLsamples.fisher$s3[f.som])
# add the metadata
mcols(s.com) = f_step3(oGRLsamples.fisher$s2[f.som], oGRLsamples.fisher$s3[f.som])

## male adult
f.male = f_step1(oGRLsamples.fisher$m2, oGRLsamples.fisher$m3)
m.com = f_step2(oGRLsamples.fisher$m2[f.male], oGRLsamples.fisher$m3[f.male])
mcols(m.com) = f_step3(oGRLsamples.fisher$m2[f.male], oGRLsamples.fisher$m3[f.male])

## female adult
f.female = f_step1(oGRLsamples.fisher$f2, oGRLsamples.fisher$f4)
f.com = f_step2(oGRLsamples.fisher$f2[f.female], oGRLsamples.fisher$f4[f.female])
mcols(f.com) = f_step3(oGRLsamples.fisher$f2[f.female], oGRLsamples.fisher$f4[f.female])

oGRLsamples.fisher.pooled = GRangesList(s=s.com, f=f.com, m=m.com)

md = list(desc='Pooled peaks common in each lifecycle stage with fisher p.values', date=paste(date()))
metadata(oGRLsamples.fisher.pooled) = md
n = make.names(paste('oGRLsamples.fisher.pooled', date(), '.rds'))
dir.create('Objects', showWarnings = F)
save(oGRLsamples.fisher.pooled, file=paste('Objects/', n, sep=''))

## create a gff for each sample
for (i in 1:length(oGRLsamples.fisher.pooled)){
  r = oGRLsamples.fisher.pooled[[i]]
  mc = as.data.frame(mcols(r))
  mc = data.frame(apply(mc, 2, function(x) signif(x, 3)))
  mcols(r) = mc
  n = make.names(paste('peaks', names(oGRLsamples.fisher.pooled)[i], date(), 'gff'))
  dir.create('Results_scripts/Gff_files', showWarnings = F)
  export(r, con=paste('Results_scripts/Gff_files/', n, sep=''), format = 'gff3')
  n = make.names(paste('peaks', names(oGRLsamples.fisher.pooled)[i], date(), 'csv'))
  write.csv(as.data.frame(r), file=paste('Results_scripts/', n, sep=''))  
}







