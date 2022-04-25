## RiboStream counting and normalisation
# This script performs a default reads counting and normalisationfor Ribo-Seq and
# RNA-Seq reads processed with RiboStream pipeline

# Load packages
library(tidyverse, quietly=T, warn.conflicts=F)
library(GenomicFeatures, quietly=T, warn.conflicts=F)
library(GenomicRanges, quietly=T, warn.conflicts=F)
library(GenomicAlignments, quietly=T, warn.conflicts=F)
library(rtracklayer, quietly=T, warn.conflicts=F)

# Load parameters and metadata
DIR = getwd()
CONFIG = "config.ini"
GTF = gsub("[']", "", gsub("GTF=", "", system(paste("cat", CONFIG, " | grep 'GTF' "), intern = TRUE)))
RIBO_LEN_MIN=as.numeric(gsub("RIBO_LEN_MIN=", "", system(paste("cat", CONFIG, " | grep 'RIBO_LEN_MIN=' "), intern = TRUE)))
RIBO_LEN_MAX=as.numeric(gsub("RIBO_LEN_MAX=", "", system(paste("cat", CONFIG, " | grep 'RIBO_LEN_MAX=' "), intern = TRUE)))

meta = read_delim("sampleSheet_alignment.csv", delim = ",")

# Helper functions
## Function to collapse a list of transcript by gene IDs
byGene <- function(GRangesList, txdb, reduce = T){
  cleaned <- unlist(GRangesList)

  # Tx to Gene ID mapping
  txNames <- names(cleaned)
  mapID <- suppressMessages(AnnotationDbi::select(txdb,
                                                  keys = as.character(txNames),
                                                  keytype = 'TXNAME',
                                                  columns = 'GENEID'))
  # Match Tx Id to Gene ID
  toKeep <- mapID$GENEID[!is.na(mapID$GENEID)]
  cleaned <- split(cleaned[!is.na(mapID$GENEID)], toKeep)

  if (reduce){
    cleaned <- IRanges::reduce(cleaned)
  }

  return(cleaned)
}

# Vectorised splitvec
splitvec <- function(vector, split, select, merge = "_"){
  processed <- sapply(vector, function(x){
    separated <- unlist(strsplit(x, split = split))[select]
    if (length(separated) > 1){
      return(paste(separated, collapse = merge))
    } else
      return(separated)
    })
  processed <- unname(processed)
  return(processed)
}

# Wrapper to organise count table
makeTab <- function(overlaps, class){
  tab <- data.frame(name = names(overlaps),
                    class = as.numeric(overlaps),
                    stringsAsFactors = F)
  colnames(tab) <- c("name", class)
  return(tab)
}


## Analysis
# Create transcript models
gtfTab <- rtracklayer::import(GTF, "GTF") %>%
  as.data.frame() %>%
  filter(transcript_type == "protein_coding"
          & tag == "CCDS")

txdb <- makeTxDbFromGFF(GTF, format = "gtf")
exons <- exonsBy(txdb, by = "gene")
cds <- cdsBy(txdb, by = "tx", use.names = T)
cds <- cds[names(cds) %in% gtfTab$transcript_id]
cdsTrimmed <- pmapFromTranscripts(IRanges(start = 30,
                                       width = sum(width(cds))-60), cds)
cdsTrimmed <- byGene(cdsTrimmed, txdb, reduce = T)
cdsFull <- byGene(cds, txdb, reduce = T)

# Create a list of bam files
bams <- list.files("bam", "*unique.bam$", full.names = T)
bamTab <- data.frame(bam = bams,
                           stringsAsFactors = F) %>%
  mutate(Experiment = case_when(grepl("prealigned", bam) ~ "Ribo",
                                TRUE ~ "RNA"),
        fileName = case_when(Experiment == "Ribo" ~ gsub("_prealigned_unique.bam", ".fastq.gz", bam),
                             TRUE ~ gsub("_unique.bam", ".fastq.gz", bam)),
        fileName = gsub("bam/", "", fileName)) %>%
  left_join(meta, by = "fileName")

# Iterating through bam files
CDScountsFull <- data.frame()
CDScountsTrim <- data.frame()
exonCounts <- data.frame()

for (b in 1:length(bamTab$bam)){
  sampleID <- gsub("_unique.bam", "", splitvec(bamTab$bam[b], "/", 2))
  samplesInfo <- bamTab %>% dplyr::filter(bam == bamTab$bam[b])

  print(paste("Processing", sampleID, ":", b, "/", length(bamTab$bam)))
  print(bamTab$bam[b])
  bf <- readGAlignments(bamTab$bam[b])
  print("Here1")

  if (samplesInfo$Experiment[1] == "Ribo"){
     bf <- bf[qwidth(bf) %in% c(28:30)]
     bf <-  GenomicAlignments::qnarrow(bf, start = 12, width = 1)
     print("Here2")
  } else {
     bf <- bf[qwidth(bf) == 50]
     bf <-  GenomicAlignments::qnarrow(bf, start = 25, width = 3)
     print("Here2")

  }

  countsTabFull <-  countOverlaps(cdsFull,  bf, minoverlap = 1, ignore.strand = T) %>%
  makeTab(class = sampleID)
  print("Here3")

  countsTabTrim <-  countOverlaps(cdsTrimmed,  bf, minoverlap = 1, ignore.strand = T) %>%
  makeTab(class = sampleID)
  print("Here4")

  exonTab <-  countOverlaps(exons,  bf, minoverlap = 1, ignore.strand = T) %>%
  makeTab(class = sampleID)
  print("Here5")

 if (b == 1){
   CDScountsFull <- countsTabFull
   CDScountsTrim <- countsTabTrim
   exonCounts <- exonTab
 } else {
   CDScountsFull <- left_join(CDScountsFull, countsTabFull, by = "name")
   CDScountsTrim <- left_join(CDScountsTrim, countsTabTrim, by = "name")
   exonCounts  <- left_join(exonCounts, exonTab, by ="name")
 }
}
print("Check1")

# Saving the output per experiment
conds <- unique(bamTab$condition)
for (c in 1:length(conds)){
  select <- which(bamTab$condition == conds[c])
print(select)
  # Raw counts
  exon_toSave <- as.data.frame(exonCounts[,c(1, select)])
  CDSfull_toSave <- as.data.frame(CDScountsFull[,c(1, select)])
  CDStrim_toSave <- as.data.frame(CDScountsTrim[,c(1, select)])
print("Check2")
  write_csv(exon_toSave, paste0("counts/", conds[c], "_Exon_counts.csv"))
  write_csv(CDSfull_toSave, paste0("counts/",conds[c], "_fullCDS_counts.csv"))
  write_csv(CDStrim_toSave, paste0("counts/",conds[c], "_trimCDS_counts.csv"))

}
print("Check3")
# Saving counting stats
countsStats <- tibble(colNames = colnames(exonCounts[-1]),
                      exonCounted = colSums(exonCounts[,-1]),
                      CDStrimCounted = colSums(CDScountsTrim[,-1]),
                      CDSfullCounted = colSums(CDScountsTrim[,-1]))

write_csv(countsStats, paste0("counts/", "Counting_stat.csv"))
