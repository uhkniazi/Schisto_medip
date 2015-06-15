# Schisto_medip
MEDIP-Seq data analysis from Schistosoma mansoni S, M and F lifecycle stages

# 00_header.R
header file to set libraries, source files and set global variables and paths

# 01_dismiss_to_gff.R
one of the first scripts in the series. used to generate an initial profile of the medip data after being processed by dismiss. 
it will generate some plots based on the peak distribution (ie. counts) in the three lifecycle stages. furthermore, based on the
initial filtering criteria, it will choose peaks (if they are present in 2 of 2, or 1 of 2 replicates) and further
subdivide them into groups i.e. peaks present in 2 lifecycle stages 1 etc. it also generates a gff which can be used
for further analysis in a genome browser and an object file used for next scripts.

# Schisto_medip_classes.R
contains functions and classes used by the script.
CFeatures - a class inherited from GRanges, and holds the feature information for the genome i.e. exon, intron etc. this
information was extracted from the gff file from sanger Schistosoma_mansoni_v5.2.gff. 
NOTE: This has been changed a little and now using the long GFF version downloaded from ftp://ftp.sanger.ac.uk/pub/project/pathogens/Schistosoma/mansoni/genome/GFF/Smansoni_gff_21032012.tar.gz

# create_features_from_gff.R
Using the gff file from ftp://ftp.sanger.ac.uk/pub/project/pathogens/Schistosoma/mansoni/genome/GFF/Smansoni_gff_21032012.tar.gz the
features i.e. Gene and CDS (exon) are extracted. The script is only using the full chromosomes (at the moment and not the smaller bits or partial unplaced features) and reads each gff file for the corresponding chromosome as a tab separated file, extracts only the type attribute Gene or CDS. We proceed by extracting the parent ID for each exon (which is something like ID=Smp_160500.1). Create a GRanges object of the exons and split this into a GRangesList using the factor parent ID. We remove about 50ish genes that have very
large number of exons (based on the distribution of number of exons in each gene) and correspondingly very long base pair region
covered by these exons. We also remove genes with only one exon, as they do not fit a casette pattern. We create the introns by
using setdiff function, create a first exon, other exons, 2k upstream and 2k downstream regions - all in the form of GRangesList
objects. The data is saved as a list of features for use by other scripts.  
  
The second part of the script uses the gff imported earlier, to create reduced exons, genes, introns, 2k upstream, 2k downstream and
the repeats file RepBasePerpignanSma52.gff. We only use the data from the full chromosomes, and the data is saved as a lFeatures 
object to be used by later scripts.

# 02_features_overlaps_casette.R
The script uses a casette feature object created earlier (create_features_from_gff.R), Pooled Medip Peaks object created earlier 
(01_dismiss_to_gff.R) and removes peaks not belonging to any class i.e. 'none'. We create a matrix with each column representing a 
feature from the casette, the rows representing each gene from that feature, and count how many times they overlap the medip GRanges
object. The summary is reported in a bar plot and a csv file.

# 03_bis_seq_vs_medip.R  
NOTE: Changes made to this SCRIPT - Still in progress and will be documented when finished  
Data published by frank lyko downloaded as bed format. It was imported as a csv file and converted to a GRanges object. The Sequence
version 5.2 was imported as a Biostrings object (DNAStringSet) and the sequence of each Bs-Seq GRanges object extracted using the
custom function f_DNAStringSet_GRangesSequenceFromDNAStringSet. Those residues that were C were marked as + stranded ranges while 
those as G were - stranded ranges. The object was saved as GRanges object and imported directly for future use. The MeDIP data was 
imported (created earlier using 01_dismiss_to_gff.R) and peaks that should occur in the male are included in the analysis. The BS-Seq
data is broken down into quantiles of proportions (how many times it was seen / total number of times it was seen). The features 
object created earlier is used (create_features_from_gff.R) for calculating feature overlaps. The main plots of interest are 
proportion of MeDIP peak and BS-Seq 5mC distribution over the features, modelled as a multinomial distribution, and confidence 
intervals calculated by resampling from a multinomial distribution.  
  
The second part of the script extends the sizes of the BS-Seq ranges to 2, and loads the dinucleotide sequence from the Biostrings
object created earlier. The sequences from the - (minus) strand are reverse complimented and the dinucleotide frequencies are
calculated for the BS-Seq data. Similarly the frequencies are calculated for the MeDIP data. The proportions are calculated as a 
multinomial distribution, with Confidence intervals calculated via resampling.

# 03_medip_overlap_genome.R
The script starts with loading the pooled medip peaks object created earlier using 01_dismiss_to_gff.R, and the features list created
earlier using the second half of the script create_features_from_gff.R. The overlaps of features (query) vs methylation data (subject) is done and the conditional probability i.e. P(y | Theta) - where y is data and Theta is the parameter - is calculated as a 
multinomial distribution. The confidence intervals are calculated by simulation. We summarize the data as bar plots with error bars for 95% confidence interval. The medip data is also split into groups (7 groups) and distributions calculated for each one 
separately as well.


