# Name: Fanny_downstream_sequences.R
# Auth: u.niazi@imperial.ac.uk
# Date: 25/06/2015
# Desc: Collect downstream (2000 bp) sequence of every smp and save in FASTA file


source('00_header.R')
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

print(show(TxDb.Hsapiens.UCSC.hg38.knownGene))

oGRgene = genes(TxDb.Hsapiens.UCSC.hg38.knownGene, columns=c('GENEID'))
# get the downstream regions
# shift 2000 bp downstream to get downstream sequences
oGRdownstream = flank(oGRgene, width = 2000, start = F)

# how many downstream regions overlap with next gene
f = overlapsAny(oGRdownstream, oGRgene)
table(f)
oGRdownstream = oGRdownstream[!f]

# get the sequences from the bsgenome database
seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, oGRdownstream)
Biostrings::writeXStringSet(seq, 'Results_scripts/hg38_downstream.fasta', format='fasta')

m = dinucleotideFrequency(seq, as.prob = F)
m = colMeans(m)
print(round(m/sum(m), 3))

m = alphabetFrequency(seq, as.prob = F)[,c('A', 'C', 'G', 'T')]
m = colMeans(m)
print(round(m/sum(m), 3))


