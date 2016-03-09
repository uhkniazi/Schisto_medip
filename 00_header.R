# File: 00_header.R
# Auth: u.niazi@imperial.ac.uk
# DESC: header file for setting global variables and loading libraries
# Date: 23/04/2015


##### header files
if (!require(GenomicFeatures)) stop('Bioconductor library GenomicFeatures required')
if (!require(rtracklayer)) stop('Bioconductor library rtracklayer required')

source('~/Dropbox/Home/Data/R/My_Libraries/NGS_functions.R')
source('Schisto_medip_classes.R')

##### PATHS
pBam = "Data_external/Bam_files"
pPeaks = 'Data_external/Dismiss_output'
pWD = getwd()

##### global variables
gcvChromosomes = c('Schisto_mansoni.Chr_ZW', 'Schisto_mansoni.Chr_1', 'Schisto_mansoni.Chr_2', 
                   'Schisto_mansoni.Chr_4', 'Schisto_mansoni.Chr_3', 'Schisto_mansoni.Chr_6', 
                   'Schisto_mansoni.Chr_7', 'Schisto_mansoni.Chr_5')

##### global functions
## get p.value from simulated data
get.prob = function(x, f){
  # get the baseline
  cBaseline = levels(f)[1]
  ivBaseline = x[which(f == cBaseline)]
  # remaining levels/groups to compare against
  cvLevels = levels(f)[-1]
  pret = sapply(cvLevels, function(l){
    x.l = x[which(f == l)]
    x.l.m = mean(x.l)
    # calculate two sided p-value
    return(min(c(sum(ivBaseline <= x.l.m)/length(ivBaseline), sum(ivBaseline >= x.l.m)/length(ivBaseline))) * 2)
  })
  return(pret)
}