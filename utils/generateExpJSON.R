params <- commandArgs(TRUE)
sampleSheet <- read.table(paste(getwd(), params[1], sep="/"), header = TRUE, sep = ",", stringsAsFactors=FALSE)
postfix <- params[2]
configFolder <- "configs/"
dir.create(file.path(getwd()), showWarnings=FALSE, recursive=TRUE)

library("dplyr", quietly=T, warn.conflicts=F)
library("jsonlite", quietly=T, warn.conflicts=F)

# Write .json files for each experiment

JSONdata <- data.frame(SampleID=character(),
		FileName=character(),
		Experiment=character(),
		Condition=character(),
                Replicate=character(),
                PairedEnd=character(),
		Merged=character())

expIDs <- unique(sampleSheet$experiment)

## Fill data related to each experiment
for(i in 1:length(expIDs)){
  expMembers <- filter(sampleSheet, experiment == expIDs[i])
  
  newExp <- data.frame(SampleID=paste(expMembers$ID[], collapse=";"))
  newExp$FileName <- paste(expMembers$fileName[], collapse=";")
  newExp$Experiment <- expIDs[i]
  newExp$Condition <- paste(expMembers$condition[], collapse=";")
  newExp$Replicate <- paste(expMembers$replicate[], collapse=";")
  newExp$PairedEnd <- paste(expMembers$pairedEnd[], collapse=";")
  newExp$Merged <- paste(expMembers$merge[], collapse=":")  
  JSONdata <- rbind(JSONdata, newExp)
}


## Save .json files for each sample
for(i in 1:dim(JSONdata)[1]){
  fileName <- paste(configFolder, as.character(expIDs[i]), postfix, ".json", sep="")
  json <- toJSON(JSONdata[i,], pretty=TRUE)
  con <- file(fileName[1])
  writeLines(json, con)
  close(con)
}

