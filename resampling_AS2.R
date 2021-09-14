#Script calculates AS2 across human samples to set the AS2 baseline of a microbe in situ
#Approach repeatedly calculates AS2, using the majority of human transcriptomes as the gold-standard human datasets, and using a few random transcriptomes as the comparison model to see how similar the datasets are to each other

wd = "C:/Users/Gina/oralMT/Analysis/Pg/"
setwd(wd)
inputversion <- "input1" #I set this to keep track of which version of the analysis I am running; I use this to pull the correct input and it is included in all of the output file names
outputversion <- "run1" #set to whatever you would like

library(tidyverse)

#input data. Want file where columns are samples and rows are genes.
#I input data that I have already normalized using vst or rlog. Otherwise, you can just load the raw data and there is code to normalize below. You need to make sure that you always normalize all of your data together (even if you are only analyzing a small subset) because the normalization varies depending on the entire dataset
df1 <- read.table(paste(inputversion, "rlogd.dds2.csv", sep = "."), sep = ",", header = TRUE)
names(df1)[1] <- "locus_tag"
rownames(df1) <- NULL

#input metadata (list of sample names and associated metadata)
metadata_file <- read_csv(file = paste(inputversion, "SunburstMetadata.csv", sep = "."))

#define sample lists
human_list <- metadata_file %>% filter(condition == "human") %>% .$filename %>% str_replace_all("-", "_")
insitu_sample_list <- human_list

all_tested_samples_list <- names(df1)[-1]
conditions_list <-  all_tested_samples_list

#load raw counts file here. Want file where columns are samples and rows are genes.
#final_raw <- read.table("input.txt", sep = "\t", header = TRUE)
#names(final_raw)[1] <- "locus_tag"
#will then need to noramlize using rlog or VST in DESeq2

#run multiple iterations of AS2 using X samples as the gold-standard and Y samples as the "model." For example, this is currently written to run 2 as the test (and the remaining samples as the gold standard) across 100 iterations.
#set iterations using c(1:XXX)
#set number of samples to "leave out" to use as the model in each iteration using model_selfvalidation_leaveout(self_validation, Y)
#define functions below
df1_double <- c(1:100) %>% map(function(x) {dummy1 <- model_selfvalidation_leaveout(self_validation, 2)
dummy1$penalty <- abs(dummy1$penalty)
names(dummy1)[2] <- paste0("penalty_", x)
return(dummy1)})  %>% purrr::reduce(left_join)

write.table(df1_double, file = paste(outputversion, "AS2_resampled_100_leaveout2.txt", sep = "_"), sep ="\t")
#output is penalties (z-scores). In R or excel, calculate AS2 for each iteration then average iterations

#################################################################################################################
#functions to calculate AS2
score_target_vs_model <- function(Target_samplenames, Model_samplenames)
{
  DF_target <- df1  %>% select(locus_tag, Target_samplenames) %>% gather(key=sample_name, value=expression, -locus_tag) %>% group_by(locus_tag) %>% dplyr::summarize(target_mean = mean(expression), target_SD = sd(expression))
  
  df_allzscores <- DF_target %>% inner_join(df1 %>% select(locus_tag, Model_samplenames)) %>% mutate_at(.vars=vars(-locus_tag, -target_mean, -target_SD), .funs=funs((.-target_mean)/target_SD))
  
  median_modelZscoreDF <- df_allzscores %>% select(locus_tag, Model_samplenames)  %>% transmute(locus_tag, penalty = pmap_dbl(.[c(-1)], function(...)  median(c(...))))
  # below should get rid of NAs
  
  if(length(which(is.na(median_modelZscoreDF$penalty)))>0)
  {
    median_modelZscoreDF[which(is.na(median_modelZscoreDF$penalty)),]$penalty <- 0
  }  
  final_DF <- median_modelZscoreDF  
  
  return(final_DF)
}

model_selfvalidation_leaveout <- function(sample_DF1, left_out_number)
{
  test_sputum <- sample(1:length(insitu_sample_list), left_out_number, replace = FALSE)
  train_sputum <- setdiff(1:length(insitu_sample_list), test_sputum)
  
  dummy1 <- score_target_vs_model(insitu_sample_list[train_sputum], insitu_sample_list[test_sputum])
  return(dummy1)
}
