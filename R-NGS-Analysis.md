---
title: "R-NGS-Analysis"
author: "Jana Gagacheva"
output:
  pdf_document: default
  html_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

The purpose of this entry is to showcase the pipeline for analysis of NGS data. Namely, the goal is to identify the guide RNAs present in each sample. Additionally, I will be determining whether there is any potential cross-contamination between samples, as this can influence the phenotypic outcome in a screen. 

Part I - handling the raw data and getting it to a point where I the sequences are trimmed and ready to be aligned to a reference file. This portion of the code is primarily done in Bash.

The first step is to demultiplex the raw data via the barcodes used in the sequencing primers. For demultiplexing, I use bcl2fastq. 

```
ulimit -n 4000

bcl2fastq -R 210409_MN00141_0035_A000H3GGFC -p 36 --output-dir demultiplexed --no-lane-splitting --sample-sheet Barcodes.csv --barcode-mismatches 1 --minimum-trimmed-read-length 20 --mask-short-adapter-reads 0
```

After the files are demultiplexed, the next step is to trim the sequences to the desired length. I use bbduk.sh from the bbmap to perform the alignment and trimming. Since I usually have a file for each sample and have many samples at once, I trim all the files using a loop in bash. In the loop, I look for all the files that end with "fastq.gz" to identify the right files. After the trimming is done, I make a new copy of each trimmed file that is re-named and ends with "_trmd.fastq.gz". In this case, I trimmed starting at base two all the way to base 19. 

```
files=/Users/jana/miniseq-run/demultiplexed
for file in $files/*.fastq.gz; do
   /Users/jana/bbmap/bbduk.sh -Xmx5g in=$file out=${file/_R1_001.fastq.gz/_trmd.fastq.gz} ftl=2 ftr=19 
done
```
The next step is to utilize FASTQC, which gives an informative output of the quality of the sequences. This step is important to gauge the quality of the samples and sequencing run. The quality control helps me become aware if I need to handle low diversity or any other sequencing issues in my sequenced sample. Again, I handle this process via a loop, where I create a QC output for every fastq file. 

```
#take a look at the quality of the trimmed sequences. 
mkdir /Users/jana/miniseq/demultiplexed/FASTQC
for file in *_R1_001.fastq.gz; do /Users/jana/Downloads/FastQC/fastqc "$file" --outdir=/Users/jana/miniseq/demultiplexed/FASTQC --extract; done
```

To be able to read all the files as sequences in R, I make a copy of each fastq.gz file to a text file. 
```
#change all the files to txt in command line
mkdir /Users/jana/demultiplexed/txts
gunzip *_trmd.fastq.gz
for file in *_trmd.fastq; do cp -v "$file" txts/"${file/_trmd.fastq/.txt}"; done
```


Part II - align the sequences to a reference file and determine the quality of each sample. This portion of the code is done in R.
```{r}
library(Biostrings)
library(ShortRead)
library(tidyverse)
library(data.table)
library(rlist)
library(forestmangr)
library(purrr)
library(stringr)
library(dplyr)
```

I first load a table that is a reference file that contains all my guide RNAs and where they are supposed to be located in my samples. 
```{r}
ref_table <- read.csv("/Users/jana/new_master.csv", header = TRUE)
original_table <- read.csv("/Users/jana/new_master.csv", header = TRUE)

#change the ref table to only contain the guides in sample #3
ref_table  = ref_table %>% filter(str_detect(Gene.Coordinates, "HG_sample3")) 

#remove the first two letters from the guides and store the result in a table called "ref_table"
##this is because i only sequenced one sample and many sequences are staring w/ a G, which reduced diversity and quality. For this reason I trimmed the sequences shorter. 
ref_table = ref_table %>% mutate(Sequence = substr(Sequence, 3, 20))
View(ref_table)
```

In the loop below, I read all the txt files from each sequence in R. I create a two-column table that contains the sample name in one column and the corresponding sequence in the other. 
For each sample, I find the guide RNA that is most common and calculate the percentage occurrence. After that, I look if there is a second guide RNA that is present in each sample. At the end, I generate a table with the guide RNAs that are in my reference file, but not present in my sequenced samples. 
```{r}
#Function to get the top 5 most common strings:
freqfunc <- function(x, n){
  tail(sort(table(unlist(x[[1]], use.names = FALSE))), n)
}

#ref_table <- tst
result_list <- list()
second_list <- list()

#############################################################################
directory <- ("/Users/jana/miniseq/demultiplexed/txts")
usingiMac <- TRUE
file_names <- list.files(directory)
#file_names <- "B3_D11_S15.txt" #file_names[[2]]
for (file_name in file_names) {
  table_file <- file.path(directory, file_name) #Make a table w/ all the file names
  
  file_name <- tools::file_path_sans_ext(basename(table_file))
  
  a1a5_txt <- read.table(table_file, sep="\t")
  a1a5 <- data.frame(a1a5_txt[lag(a1a5_txt$V1,1) %like% "@MN00141",])
  colnames(a1a5) <- c("Sequence")
  
  mostCommon_A1_A5 <- data.frame(unclass(freqfunc(a1a5, 1)))
  if (usingiMac) {
    mostCommon_A1_A5 <- cbind(Sequence = rownames(mostCommon_A1_A5), mostCommon_A1_A5)
    rownames(mostCommon_A1_A5) <- 1:nrow(mostCommon_A1_A5)
    colnames(mostCommon_A1_A5) <- c("Sequence", "Freq")
  }
  
  mostCommon_A1_A5[,"Total Number of Reads"] <- NA
  for (i in 1:nrow(mostCommon_A1_A5)) {
    mostCommon_A1_A5[[i, "Total Number of Reads"]] <- paste(nrow(a1a5))
  }
  
  random <- anti_join(a1a5, ref_table, by="Sequence")
  mostCommon_A1_A5[,"Non-library reads"] <- NA
  for (i in 1:nrow(mostCommon_A1_A5)) {
    mostCommon_A1_A5[[i, "Non-library reads"]] <- paste(nrow(random)) 
  }
  
  mostCommon_A1_A5[,"Most Common Sequence MiniSeq"] <- NA #create a column with name "Name"
  for (i in 1:nrow(mostCommon_A1_A5)) {
    mostCommon_A1_A5[[i, "Most Common Sequence MiniSeq"]] <- paste(file_name) 
  }
    
  mostCommon_A1_A5[,"Most Common Percent Occurance"] <- NA #create a column with name "Name"
  if (toString(mostCommon_A1_A5$Freq) == "") {
    print(file_name)
    print("is a degenerate")
  }
  if (nrow(a1a5) == 0 | toString(mostCommon_A1_A5$Freq) == "") {
    print("it's 0")
      percent <- 0
      for (i in 1:nrow(mostCommon_A1_A5)) {
        mostCommon_A1_A5[[i, "Most Common Percent Occurance"]] <- paste(percent) #take the row number and add it to the name of each sample in the new column "Name"
      }
  }
  else {
    percent <- c((mostCommon_A1_A5$Freq)/nrow(a1a5)*100)

    for (i in 1:nrow(mostCommon_A1_A5)) {
      mostCommon_A1_A5[[i, "Most Common Percent Occurance"]] <- paste(percent) #take the row number and add it to the name of each sample in the new column "Name"
    }
  }
  
  colnames(mostCommon_A1_A5) <- c("Sequence", "Most Common Freq", "Total Number of Reads","Non-library reads", "MiniSeq Well&PlateBC")
  
  #   mostCommon_A1_A5 <- data.frame(unclass(freqfunc(a1a5, 1)))
  #   if (usingiMac) {
  #   mostCommon_A1_A5 <- cbind(Sequence = rownames(mostCommon_A1_A5), mostCommon_A1_A5)
  #   rownames(mostCommon_A1_A5) <- 1:nrow(mostCommon_A1_A5)
  #   colnames(mostCommon_A1_A5) <- c("Sequence", "Freq")
  # }
  
  second <- data.frame(unclass(freqfunc(a1a5, 2)))
  second <- cbind(Sequence = rownames(second), second)
  rownames(second) <- 1:nrow(second)
  colnames(second) <- c("Sequence", "Freq")
  
  second_mostCom <- data.frame(second[1,])
  second_mostCom[,"Name"] <- NA #create a column with name "Name"
  for (i in 1:nrow(second_mostCom)) {
    second_mostCom[[i, "Name"]] <- paste(file_name) 
  }

  colnames(second_mostCom) <- c("Sequence", "Freq", "Name")
  second_mostCom <- merge(second_mostCom, ref_table, by="Sequence")
  second_mostCom <- second_mostCom[,c(2,3,4)]
  colnames(second_mostCom) <- c("Second Most Common Freq","MiniSeq Well&PlateBC", "Second Most Common gene")
  
  first_result <- merge(ref_table, mostCommon_A1_A5, by = "Sequence")
  result_list <- list.append(result_list, first_result) # Add our data frame to the list
  second_list <- list.append(second_list, second_mostCom) # Add our data frame to the list
}
##############################################################################
result_df = map_dfr(result_list, data.frame, stringsAsFactors=TRUE) # Convert the list back into a dataframe [as rows]
second_df <- map_dfr(second_list, data.frame, stringsAsFactors=TRUE)
result <- merge(result_df, second_df, by = "MiniSeq.Well.PlateBC", all = TRUE)
missing_sgRNA <- anti_join(ref_table, result, by="Sequence")
View(missing_sgRNA)
```

Rename the names of the barcodes and re-name the sample w/ appropriate sample nomenclature. 
```{R}
result <- data.frame(lapply(result_blue, function(x) {
  gsub("_D11_", " HG_sample3 ", x)
}))
#Remove anything following the S which is an extension given by the Miniseq machine, but I don't need.  
result <- data.frame(lapply(result_blue, function(x) {
  gsub("(.*) S.*", "\\1", x)
}))
```

Look for hits:
```{R}
#remove the rows w/ NAs in the column "Sequence"
result_filtd <- result[!is.na(result$Sequence), ]

#remove the last row called undetermined:
result_filtd = head(result_filtd, -1)

#filter only for the guides that are matching the sample in the reference file:
result_filtd_matched = result_filtd[result_filtd$Gene.Coordinates == result_filtd$MiniSeq.Well.PlateBC, ]
View(result_filtd_matched)

#transform the columns containing numbers to numeric to be able to calculate the second.most.common.percent.occurrence
result_filtd_matched <- transform(result_filtd_matched, Total.Number.of.Reads = as.numeric(Total.Number.of.Reads), 
               Second.Most.Common.Freq = as.numeric(Second.Most.Common.Freq)) 
result_filtd_matched$Second.Most.Common.Freq[is.na(result_filtd_matched$Second.Most.Common.Freq)] <- 0 #sub any value of NA to 0 in the second.most.common.frequency column
result_filtd_matched$Second.Most.Common.Percent.Occurance <- with(result_filtd_matched, (Second.Most.Common.Freq/Total.Number.of.Reads)*100 )#add a second.most.common.percent column
```


Clean up:
```{r}

#subset all the hits that have more than 85% of the most common guide (column 14)
threshhold_mostCommon = 85
result_blue_above85 = subset(result_filtd_matched, result_filtd_matched [ , 14] > threshhold_mostCommon) 
result_blue_below85 = subset(result_filtd_matched, result_filtd_matched [ , 14] < threshhold_mostCommon)

#export the "clean" results:
write.csv(result_clean, "/Users/jana/Miniseq/guides_CLEAN.csv")
```
