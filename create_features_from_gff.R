# File: create_features_from_gff.R
# Auth: u.niazi@imperial.ac.uk
# DESC: uses the gff information to create a features object
# Date: 01/05/2015


source('00_header.R')
# names of chromosomes to read data from
csChrom.names = c("Schisto_mansoni.Chr_1",  "Schisto_mansoni.Chr_2",  "Schisto_mansoni.Chr_3",
                  "Schisto_mansoni.Chr_4",  "Schisto_mansoni.Chr_5", "Schisto_mansoni.Chr_6",
                  "Schisto_mansoni.Chr_7",  "Schisto_mansoni.Chr_ZW")
csFile.names = paste('Data_external/Gff/1/', csChrom.names, '.gff.gz', sep='')
dfGff = NULL
# using the gff file provided at - see help for details
# merge the data frame to read only gene and CDS features
for (i in 1:length(csFile.names)){
  gff = read.csv(csFile.names[i], header = F, sep = '\t', skip=2)
  f = (gff$V3 == 'gene') | (gff$V3 == 'CDS')
  dfGff = rbind(dfGff, gff[f,])  
}
gc()
# remove non standard smp names
i = grepl('^ID=Smp_\\d+\\.?\\d*[:\\w]*;.+', dfGff$V9, perl=T)
dfGff = dfGff[i,]
# extract the ids for the genes and exons i.e. Smps
csID = as.character(gsub('^ID=(Smp_\\d+\\.?\\d*[:\\w]*);.+', '\\1',dfGff$V9, perl=T))
dfGff$csID = csID
# used later in the script
dfGff.bkup = dfGff

# each exon/CDS is assigned a parent gene id, so use that and make GRanges from the exons
dfGff.exon = dfGff[dfGff$V3 == 'CDS',]
# use this factor to split the GRanges object into a GRanges list based on the parent ID
fExons = gsub('(Smp_\\d+\\.?\\d*)[:\\w]*', '\\1', dfGff.exon$csID, perl=T)
# create GRanges object for exons
oGRexon = GRanges(as.character(dfGff.exon$V1), IRanges(dfGff.exon$V4, dfGff.exon$V5), strand = as.character(dfGff.exon$V7))
oGRexon$Parent = fExons 
# split it into GRangestList 
oGRLexon = split(oGRexon, fExons)
# remove genes with only one exon, large number of exons
# large total width of exons
w = sum(width(oGRLexon))

# model the cutoff
X = log(w)
r = range(X)
r[1] = floor(r[1]); r[2] = ceiling(r[2])
s = seq(r[1]-0.5, r[2]+0.5, 1)
hist(X, prob=T, main='Width of Genes', breaks=s, xlab='log Width', ylab='', ylim=c(0.0, 0.5))
s = seq(r[1], r[2], length.out = 100)
dn = dnorm(s, mean(X), sd(X))
lines(s, dn, type='l')
# choose the cutoff
c = qnorm(0.01, mean = mean(X), sd(X), lower.tail = F)
points(c, 0.0, col='red')
cutoff = which(X > c)

# remove these genes
temp = as.list(oGRLexon)
temp[names(cutoff)] = NULL
oGRLexon = GRangesList(temp)

# which genes have only one exon
X = elementLengths(oGRLexon)
cutoff = which(X == 1)
temp = as.list(oGRLexon)
temp[names(cutoff)] = NULL
oGRLexon = GRangesList(temp)

# create the gene body
oGRLgene = range(oGRLexon)
# create the introns by using setdiff, however the first intron and 1st exon may not be in the
# same order, so confirm that if we need to use this at some point
oGRLintron = vector('list', length=length(oGRLgene))
for (i in 1:length(oGRLintron)){
  oGRLintron[[i]] = setdiff(oGRLgene[[i]], oGRLexon[[i]])
}
names(oGRLintron) = names(oGRLgene)
oGRLintron = GRangesList(oGRLintron)

# create the first exon
temp = sapply(oGRLexon, function(x) x[1])
oGRLexon.1st = GRangesList(temp)
# remove first exon from other exons
temp = sapply(oGRLexon, function(x) x[-1])
oGRLexon.others = GRangesList(temp)

# create the upstream and downstream 2k regions
oGRLupstream = flank(oGRLgene, width = 2000, start = T)
oGRLdownstream = flank(oGRLgene, width = 2000, start = F)

d = paste(date(), 'gene casette created from gff file')

lFeatures = list(gene=oGRLgene, fst.exon = oGRLexon.1st,
                 exons.others=oGRLexon.others,
                 exons.all=oGRLexon,
                 introns=oGRLintron,
                 upstream=oGRLupstream,
                 downstream=oGRLdownstream,
                 desc=d)

n = paste('Objects/lGeneCasette', make.names(date()), sep='')
dir.create('Objects', showWarnings = F)
save(lFeatures, file=n)


######################################################################
###### create features object i.e. bins of features to get overlaps

# use variable created earlier
dfGff = dfGff.bkup
# each exon/CDS is assigned a parent gene id, so use that and make GRanges from the exons
dfGff.exon = dfGff[dfGff$V3 == 'CDS',]
# use this factor to split the GRanges object into a GRanges list based on the parent ID
fExons = gsub('(Smp_\\d+\\.?\\d*)[:\\w]*', '\\1', dfGff.exon$csID, perl=T)
# create GRanges object for exons
oGRexon = GRanges(as.character(dfGff.exon$V1), IRanges(dfGff.exon$V4, dfGff.exon$V5), strand = as.character(dfGff.exon$V7))
oGRexon$Parent = fExons 
# split it into GRangestList 
oGRLexon = split(oGRexon, fExons)

# create genes from exons granges list
oGRLgene = range(oGRLexon)
oGRgene = unlist(oGRLgene)

## create features
### gene scaffold for overlaps
oGRgene = reduce(oGRgene)
oGRexon = reduce(oGRexon)
oGRintron = setdiff(oGRgene, oGRexon)
oGRintron = reduce(oGRintron)
oGRprom = flank(oGRgene, width = 2000, start = T)
oGRds = flank(oGRgene, width = 2000, start = F)
# any promoters or downstream regoins at ends of chromosomes
f = width(oGRprom) < 2000
oGRprom = oGRprom[!f]
f = width(oGRds) < 2000
oGRds = oGRds[!f]

# remove promoters that run into a previous gene
f = overlapsAny(oGRprom, oGRgene)
oGRprom = oGRprom[!f]
# remove downstream regions running into the next gene
f = overlapsAny(oGRds, oGRgene)
oGRds = oGRds[!f]

# Create repeats file
gff = file.choose()
# load repeats file
gr.rep = import(gff)
# as the chromosome W is named W here and ZW in other places, rename W here to ZW
# simpler to create a new ranges object
mc = mcols(gr.rep)
sn = as.character(seqnames(gr.rep))
f = grepl('^Schisto_mansoni.Chr_W$', sn, perl=T)
sn[f] = 'Schisto_mansoni.Chr_ZW'
# create a new object
oGRrepeats = GRanges(sn, ranges(gr.rep), strand = strand(gr.rep))
mcols(oGRrepeats) = mc
oGRrepeats = oGRrepeats[seqnames(oGRrepeats) %in% gcvChromosomes]
oGRrepeats = reduce(oGRrepeats)
d = paste('compressed features i.e. binned, for schisto v5.2 created on', date())
## create features object to be used later
lFeatures = list(gene=oGRgene, exon=oGRexon,
                 intron=oGRintron,
                 prom=oGRprom,
                 downstream=oGRds,
                 repeats=oGRrepeats,
                 desc=d)

n = paste('Objects/lFeatures.', make.names(date()), sep='')
dir.create('Objects', showWarnings = F)
save(lFeatures, file=n)


################## addition in May 2016, for additional features
################## modified July 2016, with new non coding gff features
##### no part of script executed before this apart from the header import and setting global variables

lFeatures = f_LoadObject(file.choose())

# load one gff at a time
#gff = file.choose()
gff = "/home/uniazi/Data/R/Schisto_medip/Data_external/Gff/mirbase-mirRNA.gff3"
# load repeats file
gr.rep = import(gff)
# as the chromosome W is named W here and ZW in other places, rename W here to ZW
# simpler to create a new ranges object
mc = mcols(gr.rep)
sn = as.character(seqnames(gr.rep))
# make correct names i.e. change from short to long
pat = c(paste0('Chr_', 2:7), 'Chr_W')
rep = c(paste0('Schisto_mansoni.Chr_', 2:7), 'Schisto_mansoni.Chr_ZW')
for (x in seq_along(pat)) {
  f = sn %in% pat[x]
  sn[f] = rep[x]
}
# 
# f = grepl('^Schisto_mansoni.Chr_W$', sn, perl=T)
# sn[f] = 'Schisto_mansoni.Chr_ZW'
# create a new object
oGRrepeats = GRanges(sn, ranges(gr.rep), strand = strand(gr.rep))
mcols(oGRrepeats) = mc
oGRrepeats = oGRrepeats[seqnames(oGRrepeats) %in% gcvChromosomes]
oGRrepeats = reduce(oGRrepeats)
lFeatures$mirbase = oGRrepeats

## second gff
#gff = file.choose()
gff = "/home/uniazi/Data/R/Schisto_medip/Data_external/Gff/wormbase-tRNA.gff"
# load repeats file
gr.rep = import(gff)
# as the chromosome W is named W here and ZW in other places, rename W here to ZW
# simpler to create a new ranges object
mc = mcols(gr.rep)
sn = as.character(seqnames(gr.rep))
f = grepl('^Schisto_mansoni.Chr_W$', sn, perl=T)
sn[f] = 'Schisto_mansoni.Chr_ZW'
# create a new object
oGRrepeats = GRanges(sn, ranges(gr.rep), strand = strand(gr.rep))
mcols(oGRrepeats) = mc
oGRrepeats = oGRrepeats[seqnames(oGRrepeats) %in% gcvChromosomes]
oGRrepeats = reduce(oGRrepeats)
lFeatures$wormbase = oGRrepeats

d = paste('compressed features i.e. binned, for schisto v5.2 UPDATED on', date())
lFeatures$desc.2 = d

n = paste('Objects/lFeatures_reps.', make.names(date()), sep='')
dir.create('Objects', showWarnings = F)
save(lFeatures, file=n)

#################################################################################################
####### addition to create intergenic regions
# > w
# Schisto_mansoni.Chr_1 Schisto_mansoni.Chr_1.65476681 Schisto_mansoni.Chr_ZW Schisto_mansoni.Chr_ZW.59508269 
# 65476681                                               59508269 
# Schisto_mansoni.Chr_2 Schisto_mansoni.Chr_2.34464480   Schisto_mansoni.Chr_4 Schisto_mansoni.Chr_4.32115376 
# 34464480                                               32115376 
# Schisto_mansoni.Chr_3 Schisto_mansoni.Chr_3.27965507   Schisto_mansoni.Chr_6 Schisto_mansoni.Chr_6.20041393 
# 27965507                                               20041393 
# Schisto_mansoni.Chr_7 Schisto_mansoni.Chr_7.9743550    Schisto_mansoni.Chr_5 Schisto_mansoni.Chr_5.9386253 
# 9743550                                                9386253 

iSeqLens = c(59508269, 65476681, 34464480, 32115376, 27965507, 20041393, 9743550, 9386253)
names(iSeqLens) = gcvChromosomes

## create GRanges objects for these
gr = GRanges(gcvChromosomes, ranges = IRanges(1, iSeqLens))
oGRinter = gr
strand(oGRinter) = '+'
strand(gr) = '-'
oGRinter = append(oGRinter, gr)

## create intragenic from annotated features
oGRintra = c(lFeatures$gene, lFeatures$prom, lFeatures$downstream)
length(oGRintra)
oGRintra = reduce(oGRintra)
length(oGRintra)

## set difference to create non-overlapping intergenic regions
oGRinter = setdiff(oGRinter, oGRintra)
length(oGRinter)

## add these to the feature object
lFeatures$intra = oGRintra
lFeatures$inter = oGRinter

d = paste('compressed features i.e. binned, for schisto v5.2 UPDATED with intergenic features on ', date())
lFeatures$desc.3 = d

n = paste('Objects/lFeatures_reps_inter.', make.names(date()), sep='')
dir.create('Objects', showWarnings = F)
save(lFeatures, file=n)


#################################################################################
## break repeats into sub-categories
## execute the top of the script  first, i.e. header file and global variables only
gff = file.choose()
gr.rep = import(gff)
mc = mcols(gr.rep)
sn = as.character(seqnames(gr.rep))
f = grepl('^Schisto_mansoni.Chr_W$', sn, perl=T)
table(f)
sn[f] = 'Schisto_mansoni.Chr_ZW'
# create a new object
oGRrepeats = GRanges(sn, ranges(gr.rep), strand = strand(gr.rep))
mcols(oGRrepeats) = mc
oGRrepeats = oGRrepeats[seqnames(oGRrepeats) %in% gcvChromosomes]
dat = oGRrepeats$Target
fCategories = gsub('^Motif:(\\S+)".+', replacement = '\\1', x=dat, perl=T)
length(unique(fCategories))
temp = data.frame(table(fCategories))
summary(temp$Freq)
temp = temp[order(temp$Freq, decreasing = T),]
oGRrepeats$fCategories = fCategories
save(oGRrepeats, file='Objects/oGRrepeats.categories.rds')

