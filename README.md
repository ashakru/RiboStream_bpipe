# RiboStream
RiboStream is a compact Bpipe pipeline for processing and quality control of Ribo-Seq data, typically paired with corresponding bulk RNA-Seq samples. 


### Installation 

You can use GitHub to download the developer version of the pipeline or use the [lastest release]. To download the developer version:

```{bash}
gh repo clone ashakru/RiboStream
```
The RiboStream folder can be added to your `PATH`, so that you can access RiboStream executables from anywhere in your machine. 

### Dependencies

RiboStream is a Bpipe workflow, thus it reqires Bpipe to be available in PATH. RiboStream has been developed and tested for Bpipe 

RiboStream utilizes the following tools: 
- FASTQC
- Cutadapt
- Samtools
- Bowtie2
- STAR
- Mutliqc
- Deeptools

and R packages:
- tidyverse
- Ribowaltz
- GenomicRegions
- GenomicAlignments
- GenomicFeatures
- rtracklayer
- jsonlite
- tximport

#### RiboStream conda environment  

The easiest method to install all dependencies and R packages required by RiboStream is through replicating the original conda environment.

### Pre-requisites



### Config file  

RiboStream strores all parameters of the run in a single config file called `config.ini` which should be placed in the working directory. There are four  For a typical Ribo-Seq experiment it recommended to keep other parameters as default. 

### Workflows 

Full RiboStream workflow includes running two scripts: `RiboStream_preproc.sh` and `RiboStream_alignment.sh`. 

`RiboStream_preproc.sh` performs a standard quality check of sequenced reads with FASTQC, UMI-based deduplication (optional) and basic adapter trimming with Trimmomatic. FASTQC and trimming reports are collected into a single MultiQC report. 

 `RiboStream_alignment.sh` performs the core RiboStream analysis. Here, trimmed Ribo-Seq FASTQ files are first prealigned with Bowtie2 to a reference set of non-coding RNAs, such as rRNA, tRNA and snoRNA. Then, all samples are aligned to the reference genome with STAR. Obtained BAM files are sorted, indexed and uniquely mapping reads are selected. If needed, technical replicates can be merged, as specified in the sample sheet. 
 
By separating RiboStream workflow into two stages, it is possible to evaluate the quality of the data before the core analysis starts and perform extra steps (eg. additional reads trimming), if neccessary. 

### RiboStream customization 

### Support and future development

For any queries and issues with RiboStream workflow please use Issues in this GitHub repository.

If you would like to contibute to the RiboStream project facilitating standarisation and transparency of Ribo-Seq data analysis, please reach out contacting Joanna Krupka (jak75[at]cam.ac.uk)   

**Planned developments:**
1. Support for paired-end RNA-Seq analysis
2. Singularity container 
3. Automated differential translation analysis
4. Compatibility with public databases storing processed Ribo-Seq data
