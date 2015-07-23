# File: 06_bis_medip_dinucleotides.R
# Auth: u.niazi@imperial.ac.uk
# DESC: uses the bs-seq data (lyko) and medip data to calculate dinucleotide frequencies
# Date: 23/07/2015


# source header file
source('00_header.R')
library(Biostrings)

###########################################
### data loading
# load the bs-seq data object, unstranded one created earlier
oGR.bis = f_LoadObject(file.choose())
# load the MeDIP data object
oGRpooled = f_LoadObject(file.choose())
# load the biostrings object for genome version 5.2
seq = f_LoadObject(file.choose())

##########################################
### data formatting
# clean the MeDIP data to remove peaks from class 'none' and smaller scaffolds
# remove the medip peaks from class 'none'
f = which(oGRpooled$groups.lab == 'none')
oGRpooled = oGRpooled[-f]
oGRpooled = oGRpooled[seqnames(oGRpooled) %in% gcvChromosomes]

# assign strands to bs-seq object and keep only those with proportion over 0.9
f = which(oGR.bis$mcols.proportion >= 0.9)
oGR.bis = oGR.bis[f]
gc()
# rename all the W chromosome to ZW
r1 = oGR.bis[seqnames(oGR.bis) %in% gcvChromosomes]
r2 = oGR.bis[seqnames(oGR.bis) == 'Schisto_mansoni.Chr_W']
r1 = GRanges(as.character(seqnames(r1)), IRanges(start(r1), end(r1)), strand(r1), DataFrame(mcols(r1)))
r2 = GRanges('Schisto_mansoni.Chr_ZW', IRanges(start(r2), end(r2)), strand(r2), DataFrame(mcols(r2)))
oGR.bis = c(r1, r2)
s = unique(seqnames(oGR.bis))
for (i in 1:length(s)){
  # check base content for the chromosome
  base = f_DNAStringSet_GRangesSequenceFromDNAStringSet(seq, oGR.bis[seqnames(oGR.bis) == s[i]])
  print(s[i])
  print(colSums(alphabetFrequency(base)))
  # remove bases with A or T as those are not the right positions
  f = base == 'C'
  st = rep(NA, length(f))
  st[f] = '+'
  st[!f] = '-'
  f2 = which(seqnames(oGR.bis) == s[i])
  strand(oGR.bis[f2]) = st
}
rm(list = c('base', 'f', 'f2', 'r1', 'r2', 'st'))
gc()

# subset the seq object to main chromosomes
i = names(seq) %in% gcvChromosomes
seq = seq[i]

## create the random object to calculate background frequency
lf_win_denominator = function(width, bin.size){
  return(round(width/bin.size))
}

ivChrSizes = sapply(seq, length)
sam.size = length(oGRpooled)

lf_bin_chromosome = function(liChrSize, bin.size){
  return(f_bin_vector(1, end = liChrSize, bins = lf_win_denominator(liChrSize, bin.size)))
}

oGRrandom = GRanges()

# fill the object with random samples
for (i in seq_along(ivChrSizes)){
  lchr = names(ivChrSizes)[i]
  r = lf_bin_chromosome(ivChrSizes[i], 500)
  i = sample(1:nrow(r), size = sam.size, replace = F)
  r = r[i,]
  oGRrandom = append(oGRrandom, GRanges(lchr, IRanges(r$start, r$end)))
}

i = sample(seq_along(oGRrandom), size = sam.size, replace = F)
oGRrandom.sam = oGRrandom[i]

################################################################
###### dinucleotide frequency calculations
## for BS-Seq data
# resize the window size to 2
oGR.bis.w = resize(oGR.bis, 2, fix = 'start')
s = unique(seqnames(oGR.bis.w))
m = matrix(NA, nrow = length(s), ncol = 4)

for (i in 1:length(s)){
  # check base content for the chromosome
  seq.bs.plus = f_DNAStringSet_GRangesSequenceFromDNAStringSet(seq, oGR.bis.w[strand(oGR.bis.w) == '+' & seqnames(oGR.bis.w) == s[i]])
  seq.bs.minus = f_DNAStringSet_GRangesSequenceFromDNAStringSet(seq, oGR.bis.w[strand(oGR.bis.w) == '-' & seqnames(oGR.bis.w) == s[i]])
  print(s[i])
  seq.bs.minus = reverseComplement(seq.bs.minus)
  seq.bs = c(seq.bs.plus, seq.bs.minus)
  print(round(colMeans(dinucleotideFrequency(seq.bs, as.prob = T)), 3))
  mNucFreq.bis = dinucleotideFrequency(seq.bs, as.prob=F)[,c('CA', 'CC', 'CG', 'CT')]
  ivNucFreq.bis = colSums(mNucFreq.bis)
  m[i,] = ivNucFreq.bis
}
ivNucFreq.bis = round(colMeans(m), 0)
names(ivNucFreq.bis) = c('CA', 'CC', 'CG', 'CT')
rm(list = c('m', 'seq.bs.plus', 'seq.bs.minus', 'seq.bs', 'mNucFreq.bis', 'oGR.bis.w'))

### calculate for MeDIP data
# choose only peaks with a width upto the 0.95 quantile
c = quantile(width(oGRpooled), prob=0.95)
f = which(width(oGRpooled) <= c)
# make random sample of background equal to this size
i = sample(seq_along(oGRrandom), size = length(f), replace = F)
oGRrandom.sam = oGRrandom[i]

oGRpooled.sub = oGRpooled[f]
# get the sequence for plus, * and minus sides
s = unique(seqnames(oGRpooled.sub))
m = matrix(NA, nrow = length(s), ncol = 4)

for (i in 1:length(s)){
  # check base content for the chromosome
  seq.med.plus = f_DNAStringSet_GRangesSequenceFromDNAStringSet(seq, oGRpooled.sub[strand(oGRpooled.sub) == '+' & seqnames(oGRpooled.sub) == s[i]])
  seq.med.minus = f_DNAStringSet_GRangesSequenceFromDNAStringSet(seq, oGRpooled.sub[strand(oGRpooled.sub) == '-' & seqnames(oGRpooled.sub) == s[i]])
  print(s[i])
  seq.med.minus = reverseComplement(seq.med.minus)
  seq.med = c(seq.med.plus, seq.med.minus)
  print(round(colMeans(dinucleotideFrequency(seq.med, as.prob = T)), 3)[c('CA', 'CC', 'CG', 'CT')])
  mNucFreq.med = dinucleotideFrequency(seq.med, as.prob=F)[,c('CA', 'CC', 'CG', 'CT')]
  ivNucFreq.med = colMeans(mNucFreq.med)
  m[i,] = ivNucFreq.med
}

ivNucFreq.med = round(colMeans(m),0)
names(ivNucFreq.med) = c('CA', 'CC', 'CG', 'CT')
rm(list = c('m', 'seq.med.plus', 'seq.med.minus', 'seq.med', 'mNucFreq.med', 'oGRpooled.sub'))

### calculate for the random sequence
oGRpooled.sub = oGRrandom.sam
s = unique(seqnames(oGRpooled.sub))
m = matrix(NA, nrow = length(s), ncol = 4)

for (i in 1:length(s)){
  # check base content for the chromosome
  seq.med = f_DNAStringSet_GRangesSequenceFromDNAStringSet(seq, oGRpooled.sub[seqnames(oGRpooled.sub) == s[i]])
  print(s[i])
  print(round(colMeans(dinucleotideFrequency(seq.med, as.prob = T)), 3)[c('CA', 'CC', 'CG', 'CT')])
  mNucFreq = dinucleotideFrequency(seq.med, as.prob=F)[,c('CA', 'CC', 'CG', 'CT')]
  ivNucFreq = colMeans(mNucFreq)
  m[i,] = ivNucFreq
}

ivNucFreq.ran = round(colMeans(m),0)
names(ivNucFreq.ran) = c('CA', 'CC', 'CG', 'CT')
rm(list = c('m', 'seq.med', 'mNucFreq', 'oGRpooled.sub', 'ivNucFreq'))

####################################################################################
### stats using binomial model

### Medip vs background
# total trials
n = 1000
# expected probability for background
Expected.prop = ivNucFreq.ran / sum(ivNucFreq.ran)
Expected.count = rmultinom(1000, size = 1000, Expected.prop)
Expected.count = rowMedians(Expected.count)
# observed probability under experiment
Observed.prop = ivNucFreq.med / sum(ivNucFreq.med)
Observed.count = rmultinom(1000, size = 1000, Observed.prop)
Observed.count = rowMedians(Observed.count)

# expected values
df.E = data.frame(Expected.prop=round(Expected.prop, 2), Expected.count)
# observed values
df.O = data.frame(Observed.prop=round(Observed.prop, 2), Observed.count)

# binomial test
p.vals = rep(NA, length=nrow(df.E))
names(p.vals) = rownames(df.E)
for (i in seq_along(p.vals))
{
  nm = names(p.vals)[i]
  x = df.O[nm, 'Observed.count']
  p.vals[i] = binom.test(x, n, p = df.E[nm, 'Expected.prop'])$p.value
  print(paste(names(p.vals)[i], signif(p.vals[i], 3)))
}

p.adj = p.adjust(p.vals, method = 'BH')
as.data.frame(print(signif(p.adj, 3)))

### BS-seq vs background
# total trials
n = 1000

# observed probability under experiment
Observed.prop = ivNucFreq.bis / sum(ivNucFreq.bis)
Observed.count = rowMedians(rmultinom(1000, size = 1000, Observed.prop))

# observed values
df.O = data.frame(Observed.prop=round(Observed.prop, 2), Observed.count)

# binomial test
p.vals = rep(NA, length=nrow(df.E))
names(p.vals) = rownames(df.E)
for (i in seq_along(p.vals))
{
  nm = names(p.vals)[i]
  x = df.O[nm, 'Observed.count']
  p.vals[i] = binom.test(x, n, p = df.E[nm, 'Expected.prop'])$p.value
  print(paste(names(p.vals)[i], signif(p.vals[i], 3)))
}

p.adj = p.adjust(p.vals, method = 'BH')
as.data.frame(print(signif(p.adj, 3)))
