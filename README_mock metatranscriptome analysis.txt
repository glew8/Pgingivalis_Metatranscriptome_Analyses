Approach for mock (in silico) metatranscriptome construction and analysis from Lewin et al.

This analysis allows you to create a mock metatranscriptome from a set of genomes, and then determine which reads map to a (pan)genome of interest.

Note that small changes may be necessary depending on the source of your genome files used to build your mock metatranscriptomes. This was designed for genomes downloaded from the Human Oral Microbiome Database.


The overview of the steps are:
1. Construct mock metatranscriptome at different lengths using coding sequences of organisms in community. This can be done using art or a similar program.
2. Map mock metatranscriptomes to genome of interest (in my case, I mapped to the P. gingivalis ATCC 33277 genome and a pangenome using bowtie2)
3. Pull lines from the sam file that mapped using the following code:

for i in *.sam; do grep '^@' "$i" > "$i.mapped.sam"; cat "$i" | grep -v '^@' | awk -F "\t" '(and($2, 0x4) != 0x4)' | sort -u -k1,1 >> "$i.mapped.sam" ; done

4. Pull headers of mapped reads from the mapped.sam file in R using header_pull.R script. This produces a counts.txt file with the seqname of chromosome or scaffold (aka header) from each genome that mapped.
5. Make a header to organism lookup in R using lookup.R script. This produces an allheaders.csv file.
6. Compare the headers in the counts.txt file to the lookup file. I did this using lookup tables in excel. This will tell you which organisms had reads that mapped to your (pan)genome.
7. To pull the total number of sequences that mapped: wc -l *.sam.mapped.sam
8. To pull the total number of sequences corresponding to a taxon in the mock metatranscriptome, make a table with the seqname headers as lookup.linux.txt then run the following code (where XXX.fq is your full mock metatranscriptome): grep -c -f lookup.txt XXX.fq > realcounts.txt



