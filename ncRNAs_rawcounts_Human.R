################################################################################
## This file does:                                                            ##
## 1. Extract ncRNA reads from SPORTS output text files                      ##
## 2. Merge counts in one big table                                           ##
################################################################################



library(data.table)
library(dplyr) 
library(tidyr)
library(stringr)
library(tidyverse)



## Open all text files from SPORT output ##
files <- list.files(pattern="*.txt", full.names=TRUE, recursive=FALSE)

## Clean the files ##
for(x in files){
  temp <- read.table(x, header=TRUE) # load file
  temp <- temp[!grepl("piRNA", temp$Annotation ), ]
  temp$ID <- NULL
  temp$Sequence <- NULL
  temp$Match_Genome <- NULL
  temp$no_occ <- str_count(temp$Annotation, ";") ## If several names for same sequence
  temp$no_occ <- temp[, "no_occ"] + 1 ## Number occurences
  temp$WReads <- temp$Reads / temp$no_occ ## Weight the reads
  temp[, "WReads"] <- round(temp[, "WReads"], digits = 0)
  temp <- temp %>% 
    mutate(Annotation = strsplit(as.character(Annotation), ";")) %>% 
    unnest(Annotation)
  temp$ID <- paste(temp$Annotation, "-", temp$Length)
  temp <- aggregate(WReads~ID, temp, sum)    
  write.table(file=x, x=temp, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE) # write to file
  # })
}

imported_files <- lapply(files, function(x) {
  DT <- fread(x)
  new_colname <- x %>%
    basename %>%
    str_replace_all(c(".clipped_output.txt.*"="", "\\.txt$"=""))
  setnames(DT, old="WReads", new=new_colname)
  return(DT)
})

merged_data <- reduce(imported_files, merge, by = "ID", all = TRUE)

write.table( as.data.frame(merged_data), file="merge_ncRNAs.txt", row.names = FALSE )





