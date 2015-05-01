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
# extract the ids for the genes and exons i.e. Smps
csID = as.character(gsub('^ID=(Smp_\\d+\\.?\\d*[:\\w]*);.+', '\\1',dfGff$V9, perl=T))
dfGff$csID = csID
# each exon/CDS is assigned a parent gene id, so use that and make GRanges from the exons
dfGff.exon = dfGff[dfGff$V3 == 'CDS',]
# use this factor to split the GRanges object into a GRanges list based on the parent ID
fExons = gsub('(Smp_\\d+\\.?\\d*)[:\\w]*', '\\1', dfGff.exon$csID, perl=T)
# create GRanges object for exons
oGRexon = GRanges(as.character(dfGff.exon$V1), IRanges(dfGff.exon$V4, dfGff.exon$V5), strand = as.character(dfGff.exon$V7))
oGRexon$Parent = fExons 
# split it into GRangestList 
oGRLexon = split(oGRexon, fExons)



gff.file = file.choose()
oGR.gff = import(gff.file)


## get the genes out first
dfMcols = mcols(oGR.gff)

f = oGR.gff$type == 'gene'
oGR