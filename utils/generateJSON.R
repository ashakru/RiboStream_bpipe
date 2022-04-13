#R script to generate JSON files with settings for each Ribo-Seq sample
params <- commandArgs(TRUE)
sampleSheet <- read.table(paste(getwd(), params[1], sep="/"), header = TRUE, sep = ",", stringsAsFactors=FALSE)
postfix <- params[2]
configFolder <- "configs/"
dir.create(file.path(getwd()), showWarnings=FALSE, recursive=TRUE)

library("dplyr", quietly=T, warn.conflicts=F)
library("jsonlite", quietly=T, warn.conflicts=F)

##### Write .json files

## Initialize JSONdata object
JSONdata <- sampleSheet

## Delete NAs, fill "pairedEnd" column with "false" values
for(i in 1:dim(JSONdata)[1]){
  for(j in 1:dim(JSONdata)[2]){
    if(is.na(JSONdata[i, j])){
      JSONdata[i, j] <- ""
    }
    if(names(JSONdata)[j] == "pairedEnd" && JSONdata[i, j] == ""){
      JSONdata[i, j] <- "false"
    }
  }
}

## If structure of the .csv is the initial, change columns and fill them accordingly
if(length(colnames(JSONdata)) == 11){
  if (all(colnames(JSONdata) == c("ID","fileName","fileNamePE","experiment","condition",
  "replicate","merge","phredScore","adapter","adapterPE","readLength"))){
    ## Fill "merge" column
    merged <- list(0)
    JSONdata$merge[JSONdata$merge == ""] <- "None"

    for(i in 1:dim(JSONdata)[1]){
      if(JSONdata$merge[i] != "" && JSONdata$merge[i] != "None"){
        if(!JSONdata$ID[i] %in% merged[[1]]){
          mergeDF <- filter(JSONdata, experiment == JSONdata$experiment[i] & condition == JSONdata$condition[i] & replicate == JSONdata$replicate[i]
                            & merge == JSONdata$merge[i])
          if(dim(mergeDF)[1] > 1){
            mergeIDlist <- mergeDF$ID
            merged <- list(c(unlist(merged), unlist(mergeIDlist[2:length(mergeIDlist)])))
            JSONdata$merge[JSONdata$ID %in% mergeIDlist] <- ""
            JSONdata$merge[i] <- paste(mergeIDlist[], collapse=";")
          }
        }
      }
    }

    ## Delete unnecessary columns
    JSONdata$lane <- NULL
  }

}else{
  JSONdata$pairedEnd <- tolower(as.character(JSONdata$pairedEnd))
  JSONdata$phredScore <- 33
}

## Save .json files
for(i in 1:dim(JSONdata)[1]){
  fileName <- paste(configFolder, as.character(JSONdata$ID[i]), postfix, ".json", sep="")
  json <- toJSON(JSONdata[i,], pretty=TRUE)
  con <- file(fileName[1])
  writeLines(json, con)
  close(con)
}
