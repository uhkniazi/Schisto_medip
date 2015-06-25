# Name: Fanny_downstream_sequences.R
# Auth: u.niazi@imperial.ac.uk
# Date: 25/06/2015
# Desc: Collect downstream (2000 bp) sequence of every smp and save in FASTA file


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

# each gene is assigned a parent gene id, so use that and make GRanges from the genes
dfGff.gene = dfGff[dfGff$V3 == 'gene',]
# use this factor to split the GRanges object into a GRanges list based on the parent ID
fGene = dfGff.gene$csID
# create GRanges object for exons
oGRgene = GRanges(as.character(dfGff.gene$V1), IRanges(dfGff.gene$V4, dfGff.gene$V5), strand = as.character(dfGff.gene$V7))
oGRgene$ID = as.character(fGene)
# shift 2000 bp downstream to get downstream sequences
oGRdownstream = flank(oGRgene, width = 2000, start = F)

# how many downstream regions overlap with next gene
f = overlapsAny(oGRdownstream, oGRgene)
oGRdownstream = oGRdownstream[!f]
# check if any start or ends overrun the ends of the genome
st = start(oGRdownstream)
i = which(st < 0)
oGRdownstream = oGRdownstream[-i]

library(Biostrings)

seq = f_LoadObject(file.choose())

s = unique(as.character(seqnames(oGRdownstream)))
seq.ds.p = DNAStringSet()
seq.ds.m = DNAStringSet()

for (i in 1:length(s)){
  # get sequence for + and - strands of these regions
  seq.ds.plus = f_DNAStringSet_GRangesSequenceFromDNAStringSet(seq, oGRdownstream[strand(oGRdownstream) == '+' & seqnames(oGRdownstream) == s[i]])
  seq.ds.minus = f_DNAStringSet_GRangesSequenceFromDNAStringSet(seq, oGRdownstream[strand(oGRdownstream) == '-' & seqnames(oGRdownstream) == s[i]])
  print(s[i])
  # reverse compliment the - strand sequence
  seq.ds.minus = reverseComplement(seq.ds.minus)
  # assign names to the sequences based on their smp
  names(seq.ds.plus) = oGRdownstream[strand(oGRdownstream) == '+' & seqnames(oGRdownstream) == s[i]]$ID
  names(seq.ds.minus) = oGRdownstream[strand(oGRdownstream) == '-' & seqnames(oGRdownstream) == s[i]]$ID
  seq.ds.p = c(seq.ds.p, seq.ds.plus)
  seq.ds.m = c(seq.ds.m, seq.ds.minus)
}

# write 2 fasta file for plus and minus strands
Biostrings::writeXStringSet(seq.ds.p, 'Results_scripts/downstream_plus.fasta', format='fasta')
Biostrings::writeXStringSet(seq.ds.m, 'Results_scripts/downstream_minus.fasta', format='fasta')

m = dinucleotideFrequency(seq.ds.p, as.prob = F)
m = colMeans(m)
print(round(m/sum(m), 3))

m = dinucleotideFrequency(seq.ds.m, as.prob = F)
m = colMeans(m)
print(round(m/sum(m), 3))

m = alphabetFrequency(seq.ds.p, as.prob = F)[,c('A', 'C', 'G', 'T')]
m = colMeans(m)
print(round(m/sum(m), 3))

m = alphabetFrequency(seq.ds.m, as.prob = F)[,c('A', 'C', 'G', 'T')]
m = colMeans(m)
print(round(m/sum(m), 3))

