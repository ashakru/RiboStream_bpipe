params <- commandArgs(TRUE)
sampleSheet <- read.table(paste(getwd(), params[1], sep="/"), header = TRUE, sep = ",", stringsAsFactors=FALSE)
postfix <- params[2]
configFolder <- "configs/"
dir.create(file.path(getwd()), showWarnings=FALSE, recursive=TRUE)

library("dplyr", quietly=T, warn.conflicts=F)
library("jsonlite", quietly=T, warn.conflicts=F)

# Write .json files for each experiment

JSONdata <- data.frame(ID=character(),
											 FileName=character(),
											 Experiment=character(),
											 Condition=character(),
											 Replicate=character(),
											 Merged=character(),
										 	 Class=character())

sampleSheet <- filter(sampleSheet, merge != "")
expIDs <- unique(sampleSheet$ID)
conditions <- unique(sampleSheet$condition)
experiments <- unique(sampleSheet$experiment)
types <- c("genome","transcriptome")

## List all bam files
bams <- list.files("bam", ".bam$")
bams <- bams[!grepl("merged.bam$", bams)]


bam_tab <- data.frame(SampleName = character(),
		SampleID = character(),
		Experiment = character(),
		Condition = character())

## Add attributes for each bam file
for (b in 1:length(bams)){
	newbam <- data.frame(SampleName = bams[b])
	newbam$sampleID <- paste(unlist(strsplit(bams[b], "[.]"))[1:2], collapse = ".")
	newbam$Experiment <- sampleSheet[grepl(newbam$sampleID, sampleSheet$fileName),4]
	newbam$Condition <- sampleSheet[grepl(newbam$sampleID, sampleSheet$fileName),5]

	bam_tab <- rbind(bam_tab, newbam)
}

bam_tab$alignment_type <- types[1]
bam_tab$alignment_type[grepl("tx", bam_tab$SampleName)] <- types[2]


## Fill data related to each experiment
for(i in 1:length(conditions)){
	for (e in 1:length(unique(experiments))){
		for (t in 1:2){
			expMembers <- bam_tab %>% dplyr::filter(Condition == conditions[i] & Experiment == experiments[e] & alignment_type == types[t])

			if (all(expMembers$alignment_type == "transcriptome"))
				type <- "tx"
			if (all(expMembers$alignment_type == "genome"))
				type <- "gn"

			newExp <- data.frame(ID = paste(expMembers$sampleID[1], type, sep = "_"),
					     FileName = expMembers$SampleName[1],
					     Experiment = expMembers$Experiment[1],
					     Condition = expMembers$Condition[1],
					     Replicate = dim(expMembers)[1],
					     Class = expMembers$alignment_type[1])

			newExp$Merged <- c(paste(as.character(expMembers$SampleName), collapse = ";"))
			JSONdata <- rbind(JSONdata, newExp)

			## Save merged.txt
			merge_tab <- data.frame(samples = paste0(getwd(), "/bam/", expMembers$SampleName))
			filename <- paste0(getwd(), "/", paste(expMembers$sampleID[1], type, sep = "_"), "_merge.txt")
			readr::write_delim(merge_tab, filename,col_names = FALSE, delim = "\t")
		}
	}
}
JSONdata <- na.omit(JSONdata)

## Save .json files for each sample
for(i in 1:dim(JSONdata)[1]){
	name <- paste(JSONdata$Condition[i], JSONdata$Class[i], JSONdata$Experiment[i], sep = "_")
  fileName <- paste(configFolder, as.character(name), "_orf", ".json", sep="")
  json <- toJSON(JSONdata[i,], pretty=TRUE)
  con <- file(fileName[1])
  writeLines(json, con)
  close(con)
}
