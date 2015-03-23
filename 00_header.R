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

