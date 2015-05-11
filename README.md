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

# 02_features_overlaps_casette.R
The script uses a casette feature object created earlier (create_features_from_gff.R), Pooled Medip Peaks object created earlier 
(01_dismiss_to_gff.R) and removes peaks not belonging to any class i.e. 'none'. We create a matrix with each column representing a 
feature from the casette, the rows representing each gene from that feature, and count how many times they overlap the medip GRanges
object. The summary is reported in a bar plot and a csv file.





