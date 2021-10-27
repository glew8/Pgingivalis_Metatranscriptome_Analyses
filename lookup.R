#R script to pull headers from CDS files to identify genomes matching to coding sequences. Headers are the seqname of the genome's chromosome or scaffold
#input is folder with headers.txt files for all organisms included in mock metatranscriptome. This is made from the cds_from_genomic.fna files using the following command: for i in *cds_from_genomic.fna; do grep \> "$i" > headers/"$i.txt" ; done
#Output is lookup between headers and the name of the genome assembly (GCA)


wd = "/glewin/scratch/genomes/headers"
setwd(wd)
filenames <- dir(pattern = "_cds_from_genomic.fna.txt")
conditions <- sub("_cds_from_genomic.fna.txt","", filenames) #sub(pattern, replacement, x)
filenames
conditions

library(tidyverse)

all_headers <- data.frame(V1=c(0))
for (i in 1:length(conditions)) {
  headers <- read.csv(paste(paste(conditions[i], sep="/"), "cds_from_genomic.fna.txt", sep="_"), header = FALSE, sep = "_")
  headers2 <- headers %>% select(V1)
  headers2 <- headers2 %>% mutate(V1 = str_replace(V1, ">", ""))
  headers2$GCA <- conditions[i]
  headers3 <- headers2[!duplicated(headers2), ]
  all_headers <- merge(all_headers, headers3, all=T)
} 
all_headers <- tail(all_headers, n=-1)

write.csv(all_headers, file="allheaders.csv", row.names=FALSE)