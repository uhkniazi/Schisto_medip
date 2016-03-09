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

### get pairwise comparisons for simulated p.values with multiple test correction
get.pairwise.prob = function(x, fac, p.adjust.method='none'){
  f = function(lev.1, lev.2){
    sapply(seq_along(lev.1), function(i){
      ivBaseline = x[which(fac == lev.1[i])]
      x.l = x[which(fac == lev.2[i])]
      x.l.m = mean(x.l)
      return(min(c(sum(ivBaseline <= x.l.m)/length(ivBaseline), sum(ivBaseline >= x.l.m)/length(ivBaseline))) * 2)
    })}
  # levels for the factor
  lev = levels(fac)
  # perform pairwise comparison using outer product
  out = outer(lev[-length(lev)], lev[-1L], f) 
  dimnames(out) = list(paste(lev[-length(lev)], '*', sep=''), lev[-1L])
  out = t(out)
  ## adjust the p.value for multiple testing
  out[upper.tri(out)] = NA
  out[lower.tri(out, diag = TRUE)] = p.adjust(out[lower.tri(out, diag = TRUE)], method=p.adjust.method)
  return(out)
}



pairwise.table()