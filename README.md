# Pgingivalis_Metatranscriptome_Analyses

Scripts to perform:
1. Accuracy score analyses (scripts in main folder)
2. Mock metatranscriptome analyses (scripts in labeled subfolder)

################################################################
ACCURACY SCORE SCRIPTS

Scripts to analyze the transcriptional profile of Porphyromonas gingivalis in human metatranscriptomes and laboratory model transcriptomes

AS2andSunbursts.R: this is the standard script that will calculate the AS2 for a set of transcriptomes and plot the sunburst. This script will also calculate the AS2 for individual samples

resampling_AS2.R: Script calculates AS2 across human samples to set the AS2 baseline of a microbe in situ

input1.SunburstMetadata.csv: example metadata file. Can add additional columns to analyze data in different combinations.

input1.SunburstCategories.txt: example input to make sunburst graphs using TIGRFAM categories

input1.AnnotationTIGRFAM.txt: example input to make sunburst graphs using TIGRFAM categories

INSTALL the following packages:

use install.packages command in R unless otherwise noted

1. tidyverse (confirm installation of ggplot2, readr, tibble, purrr, dplyr, tidyr)
2. cowplot
3. zeallot
4. ggsunburst: see below

ggsunburst install (This is what worked for me):

>install.packages(c("devtools", "reticulate", "reshape2", "rappdirs", "backports"))
>
>library(devtools)
>
>install_github("didacs/ggsunburst")
>
>library(reticulate)
>
>py_install("six")

################################################################
MOCK METATRANSCRIPTOME SCRIPTS

Scripts to create a mock metatranscriptome from a set of genomes, and then determine which reads map to a (pan)genome of interest. This is useful for understanding the specificity of your analyses and for determining the minimum read length that will accurately map to an organism of interest.

README_mock metatranscriptome analysis.txt: steps for performing analysis, including some code

header_pull.R: script to parse mapped.sam files with the goal of identifying the genome origin of mapped reads

lookup.R: script to make a reference table to match headers from CDS files with the corresponding genome assem to identify genomes matching to coding sequences. Headers are the seqname of the genome's chromosome or scaffold
