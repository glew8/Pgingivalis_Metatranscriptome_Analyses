#R script to parse mapped.sam files with the goal of identifying the genome origin of mapped reads.
#input is sam file parsed to only include mapped reads. Parse with shell script: for i in *Pg_pan27.sam; do cat "$i" | grep -v '^@' | awk -F "\t" '(and($2, 0x4) != 0x4)' | sort -u -k1,1 >> "$i.mapped.sam" ; done
#Output is unique genome identifier (seqname of chromosome or header) and count.


wd = "/glewin/scratch/"
setwd(wd)
filenames <- dir(pattern = ".sam.mapped.sam")
conditions <- sub(".sam.mapped.sam","", filenames) #sub(pattern, replacement, x)
filenames
conditions

library(splitstackshape)

for (i in 1:length(conditions)) {
  data <- read.csv(paste(paste(conditions[i], sep="/"), "sam.mapped.sam", sep="."), sep="\t",fill=TRUE,header = F,quote = "",col.names = paste("V",seq_len(20), sep=""))
  data <- tail(data, n=-3)
  newdata <- cSplit(data,"V1", sep="_", direction="wide", fixed=T)
  counts <- table(newdata$V1_1)
  write.table(counts, paste(paste(conditions[i], sep="/"), "counts.txt", sep="."), sep = "\t", row.names=FALSE)
  rm(counts)
}