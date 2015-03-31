# File: Schisto_medip_classes.R
# Auth: u.niazi@imperial.ac.uk
# DESC: classes and functions used by scripts
# Date: 25/11/2014


library(methods)

####### class CFeatures
# Desc: an extention of the class GRanges in package GenomicFeatures. holds some additional associated data
#       with each ranges that will provide mapping to other feature objects
#       a possible use is if we have a parent set of ranges that have children features - then using
#       slots like iID as primary keys and iParent as foriegn keys can provide mapping between two sets of 
#       related features

library(GenomicFeatures)

setClass('CFeatures', slots=list(iID='numeric', csName='character', csType='character', iParent='numeric'), 
         contains='GRanges')

CFeatures = function(feature.id, feature.type, gr.feature, feature.name='none', parent.id=0){
  # error checks
  if (length(feature.id) == 0 || length(feature.id) != length(parent.id)) stop('IDs length can not be zero')
  # assign data 
  new('CFeatures', iID=feature.id, csType=feature.type, iParent=parent.id, csName=feature.name, 
      gr.feature)
} # end constructor

# accessor functions
CFeatures.getID = function(obj) obj@iID
CFeatures.getName = function(obj) obj@csName
CFeatures.getType = function(obj) obj@csType
CFeatures.getParent = function(obj) obj@iParent

CFeatures.setParent = function(obj, parent.id) obj@iParent = parent.id

########### end class
