#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("GenomicFeatures"))
suppressPackageStartupMessages(library("tximport"))
args <- commandArgs(TRUE)

#Construct sample sheet
samplesheet <- read.csv(args[1])

#Make database of transcript and gene IDs
gtf <- args[2]
txdb <- makeTxDbFromGFF(gtf, "gtf")
k <- keys(txdb, keytype = "TXNAME")

tx2gene <- select(txdb, k, "GENEID", "TXNAME")
tx2gene$TXNAME <- sapply(tx2gene$TXNAME, function(x){unlist(strsplit(x, split = "[.]"))[1]})

#Two types of salmon run (with and without --validateMappins option)

print(paste("Importing quant.sf files ..."))
	
files <- list.files(file.path("bam"), pattern = "quant.sf", recursive = T, full.names = T)
names(files) <- sapply(files,
			function(x){
				full_name <- unlist(strsplit(as.character(x),
								split = "/"))[2]
				sampleID <- paste(unlist(strsplit(full_name,
								split = "[.]"))[1:2], collapse = ".")
				return(sampleID)
				}

)

txi_gene <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = T)
txi_tx <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = T, txOut = T )

save(txi_gene, file = "RNA_transcriptome/RNAseq_tximport_gene.RData")
save(txi_tx, file = "RNA_transcriptome/RNAseq_tximport_tx.RData")
