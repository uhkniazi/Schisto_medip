# File: 07_rna_seq_aza.R
# Auth: u.niazi@imperial.ac.uk
# DESC: rna-seq data of controls vs aza treated samples
# Date: 23/07/2015


# source header file
source('00_header.R')
library(DESeq2)
library(GO.db)

mDat = as.matrix(read.csv(file.choose(), sep = ',', header=T))
mCounts = apply(mDat[,2:7], 2, as.integer)
rownames(mCounts) = mDat[,1]

f = colnames(mCounts)
f = gl(2, k = 3, labels = c('C', 'S'))
dfDesign = data.frame(condition=f, row.names = colnames(mCounts))

# call deseq2 constructor
oDseq = DESeqDataSetFromMatrix(mCounts, dfDesign, design = ~ condition)
oDseq = DESeq(oDseq)

oRes = results(oDseq, contrast = c('condition', 'S', 'C'))
plotMA(oDseq)

dfRes = as.data.frame(oRes[which(oRes$padj < 0.1),])
n = rownames(dfRes)
mCounts.norm = counts(oDseq, normalized=T)
mCounts.norm = mCounts.norm[n,]

# match the gene names with genes from casettes
dfCasetteGroups = f_LoadObject(file.choose())
dfCasetteGroups = dfCasetteGroups[dfCasetteGroups$female,]
n2 = rownames(dfCasetteGroups)
n2 = gsub('(.+)\\.\\d', replacement = '\\1', x = n2)
# remove alternative spliced versions which will be duplicated names
i = !duplicated(n2)
dfCasetteGroups = dfCasetteGroups[i,]
# replace names without .
n2 = rownames(dfCasetteGroups)
n2 = gsub('(.+)\\.\\d', replacement = '\\1', x = n2)
rownames(dfCasetteGroups) = n2
# get these genes from the results object
dfRes.sub = dfRes[n2,]
dfRes.sub = na.omit(dfRes.sub)
dfCasetteGroups = dfCasetteGroups[rownames(dfRes.sub),]

# create groups based on response direction
fGroups = ifelse(dfRes.sub$log2FoldChange > 0, 'UP', 'DOWN')
fGroups = factor(fGroups, levels = c('DOWN', 'UP'))

# fit a regression model to check for relationships

dfData = data.frame(groups=factor(dfCasetteGroups$groups.lab), fGroups)
dfData = dfCasetteGroups[,2:6]
dfData$fGroups = fGroups
library(MASS)
library(tree)
fit.1 = glm(fGroups ~ ., data=dfData, family='binomial')
fit.2 = lda(fGroups ~., data=dfData)
fit.3 = tree(fGroups ~ ., data=dfData)

## load the data from the GOTerms
dfGO = read.csv(file.choose(), sep='\t', header = F, stringsAsFactors=F)
# choose goterms from process P
# f = dfGO$V9 == 'P'
# dfGO = dfGO[f,]
f = dfGO$V2 %in% rownames(mCounts.norm)
dfGO = dfGO[f,]

dfGraph = dfGO[,c('V2', 'V5')]
mCounts = mCounts.norm[unique(dfGraph$V2),]

source('../CGraphClust/CGraphClust.R')
fGroups = gl(2, k = 3, labels = c('C', 'S'))
mCounts = log(mCounts+1)

mCounts = t(mCounts)
mCor = cor(mCounts)

hist(sample(mCor, 1000, replace = F), prob=T)

oGr = CGraphClust(dfGraph, abs(mCor), iCorCut = 0.5)

p.old = par(mar=c(6, 3, 1,1))
mSig = (plot.significant.expressions(oGr, t(mCounts), fGroups, cex.axis=0.8))$means

#cvClust = getClusterLabels(oGr)
#cvClust = unique(cvClust)

dfGo.def = AnnotationDbi::select(GO.db, rownames(mSig), columns = c('ONTOLOGY', 'DEFINITION'), keytype = 'GOID')

plot.cluster.expressions(oGr, t(mCounts), fGroups, csClustLabel = 'GO:0005215' )
plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = 'GO:0005215')
plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = 'GO:0004222')
plot.heatmap.cluster(oGr, t(mCounts), csClustLabel = 'GO:0008380')

plot.heatmap.means(oGr, t(mCounts))

pr.out = plot.components(oGr, t(mCounts), fGroups, bStabalize = F)
biplot(pr.out)

mCent = mPrintCentralitySummary(oGr)
cvGenes = unique(unlist(lGetTopVertices(oGr)))
getLargestCliques(oGr)



