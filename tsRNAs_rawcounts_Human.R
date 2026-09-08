################################################################################
## This file does:                                                            ##
## 1. Extract tsRNA reads from SPORTS output text files                       ##
## 2. Merge counts in one big table                                           ##
################################################################################

#----[ load packages ]------#
{
  library("data.table")
  library("dplyr") 
  library("tidyr")
  library("stringr")
  library("tidyverse") 
  
}



#----[ Extract tsRNA reads from SPORTS output text files ]------#
{
  files <- list.files(pattern="*output.txt", full.names=TRUE, recursive=FALSE)
  for(x in files){
    temp <- read.table(x, header=TRUE) # load file
    temp <- temp[grep("tRNA", temp$Annotation ), ]
    temp$ID <- NULL
    temp$Sequence <- NULL
    temp <- temp[temp$Length>=15 & temp$Length<=36,]
    temp$Length <- NULL
    temp$Match_Genome <- NULL
    temp$no_occ <- str_count(temp$Annotation, ";")
    temp$no_occ <- temp[, 3] + 1
    temp$WReads <- temp$Reads / temp$no_occ
    temp[, 4] <- round(temp[, 4], digits = 0)
    temp <- temp %>% 
      mutate(Annotation = strsplit(as.character(Annotation), ";")) %>% 
      unnest(Annotation)
    temp <- aggregate(WReads~Annotation, temp, sum)    # write to file
    write.table(file=x, x=temp, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  }
  
}


#----[ Merge miRNA reads from new cleaned text files ]------#
{
  imported_files <- lapply(files, function(x) {
    DT <- fread(x)
    new_colname <- x %>%
      basename %>%
      str_replace_all(c("_S[0-9].*"="", "\\.txt$"=""))
    setnames(DT, old="WReads", new=new_colname)
    return(DT)
  })
  
  merged_data <- reduce(imported_files, merge, by = "Annotation", all = TRUE)
  
  write.table( as.data.frame(merged_data), file="All_merge_tsRNAs.txt", row.names = FALSE )
  
}







