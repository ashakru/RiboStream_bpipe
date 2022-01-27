# RiboStream
RiboStream is a compact Bpipe pipeline for processing and quality control of Ribo-Seq data, typically paired with corresponding bulk RNA-Seq samples. 

### Installation 
You can use GitHub to download the developer version of the pipeline or use the [lastest release]. To download the developer version:

```{bash}
gh repo clone ashakru/RiboStream
```
The RiboStream folder can be added to your `PATH`, so that you can access RiboStream executables from anywhere in your machine. You can do it, for example, following this [tutorial](https://www.howtogeek.com/658904/how-to-add-a-directory-to-your-path-in-linux/)

### Dependencies

RiboStream is a Bpipe workflow, thus it reqires Bpipe to be available in PATH. RiboStream has been developed and tested for Bpipe 

RiboStream utilizes the following tools: 
- [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [trimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
- [Samtools](https://www.htslib.org/doc/samtools.html)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [STAR](https://github.com/alexdobin/STAR)
- [MutliQC](https://multiqc.info)
- [Deeptools](https://deeptools.readthedocs.io/en/develop/)

and R packages:
- [tidyverse](https://www.tidyverse.org)
- [Ribowaltz](https://github.com/LabTranslationalArchitectomics/riboWaltz)
- [GenomicAlignments](https://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html)
- [GenomicFeatures](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
- [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html)
- [jsonlite](https://cran.r-project.org/web/packages/jsonlite/index.html)
- [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html)
- [Salmon](https://combine-lab.github.io/salmon/getting_started/)

#### RiboStream conda environment  

The easiest method to install all dependencies and R packages required by RiboStream is through replicating the original conda environment. `RiboStream.yaml` file with all the dependencies and versions is available in 'env' folder. To install conda follow [this tutorial](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) and to recreate the RiboStream environment simply:

```{bash}
conda env create -f RiboStream.yml
```

### Pre-requisites

#### 1. Input FASTQ files 
Demultiplexed FASTQ files should be placed in a `fastq` folder in the working directory. Alternatively, a list of SRA accession numbers can be provided in the sample sheet, see below. 

#### 2. Sample sheet with samples metadata
A sample sheet is a `.csv` file containing all importat information about each sample to ensure right workflow execution. An example file, `sampleSheet_test.csv` can be found in `test\` folder in this GitHub repository. The sample sheet file consists of 6 columns: 
| Column name   | Explanation                                                                                   |
| ------------- | --------------------------------------------------------------------------------------------- |
| ID            | Unique sample ID (if processing from FASTQ files) OR SRA accession number                     |
| fileName      | Name of a FASTQ files in the `fastq` folder, leave empty if SRA accession number was provided |
| experiment    | 'Ribo' for Ribo-Seq or 'RNA' for RNA-Seq                                                      |
| group         | Group name that will be used to return separate count matrices and QC plots per group.        |
| condition     | Experimental condition, eg. control or shRNA                                                  |
| replicate     | Biological replicate number                                                                   |
| merge         | Indicates if technical relicates are to be merged together, leave empty if not required       |
| adapter       | Adapter sequence passed to `-a` TrimGalore argument                                           |
| trimQuality   | Trimming quality passed to `-q` TrimGalore argument                                           |


#### 3. Bowtie2 and STAR index
Next step is to create Bowtie2 and STAR indexes for the reference genome. They can be stored anywere you like (not necessarily in your working directory where all your FASTQ files are). The location of the indexes will be specified in the config file. 

##### Bowtie2
Bowtie2 is used to prealign trimmed Ribo-Seq reads to a set of non-coding RNAs. The advantage of this extra data purification step is an increase in the proportion of true ribosome-derived mRNA footprints in the final BAM file. This facilitates the interpretation of characteristic periodical alignment pattern of Ribo-Seq reads and improves the estimation of the ribosomal P-site location. Additional benefit is smaller size of the final output file that accelerates further processing. Reference sequences of rRNA, tRNA and snoRNA can be downloaded from RefSeq, see an example FASTA file in `test/ref` folder in this GitHub repository. To create Bowtie2 index plese refer to the original Bowtie2 manual. An example command used to create Bowtie2 index for RiboStream tests was:
```{bash}
bowtie2 index 
```

##### STAR
RiboStream uses STAR to align both Ribo-Seq and RNA-Seq read to the reference genome. You can build your STAR index with any reference sequence. For example, reference genome for human and mouse can be downloaded from [Gencode](https://www.gencodegenes.org). An example command used to create STAR index for RiboStream tests was:
```{bash}
STAR index 
```

#### 5. GTF file with reference gene models  
A GTF file containing exon and CDS genomic coordinates is required by STAR and at several steps of the routine Ribo-Seq quality check. For example, RiboStream was developed and tested with a Comprehensive gene annotation (CHR) downloaded from [Gencode](https://www.gencodegenes.org/human/)

#### 4. `config.ini` file with tools parameters 
RiboStream stores all parameters of the run in a single config file called `config.ini` which should be placed in the working directory. There are four  For a typical Ribo-Seq experiment it recommended to keep other parameters as default. 

### Workflows 

Full RiboStream workflow includes running two scripts: `RiboStream_preproc.sh` and `RiboStream_alignment.sh`. 

`RiboStream_preproc.sh` performs a standard quality check of sequenced reads with FASTQC, UMI-based deduplication (optional) and basic adapter trimming with Trimmomatic. FASTQC and trimming reports are collected into a single MultiQC report. 

 `RiboStream_alignment.sh` performs the core RiboStream analysis. Here, trimmed Ribo-Seq FASTQ files are first prealigned with Bowtie2 to a reference set of non-coding RNAs, such as rRNA, tRNA and snoRNA. Then, all samples are aligned to the reference genome with STAR. Obtained BAM files are sorted, indexed and uniquely mapping reads are selected. If needed, technical replicates can be merged, as specified in the sample sheet. 
 
By separating RiboStream workflow into two stages, it is possible to evaluate the quality of the data before the core analysis starts and perform extra steps (eg. additional reads trimming), if neccessary. 

### Compatibility with HPC

### RiboStream customization 

### Support and future development

For any queries and issues with RiboStream workflow please use Issues in this GitHub repository.

If you would like to contibute to the RiboStream project facilitating standarisation and transparency of Ribo-Seq data analysis, please reach out contacting Joanna Krupka (jak75[at]cam.ac.uk)   

**Planned developments:**
1. Support for paired-end RNA-Seq analysis
2. Singularity container 
3. Automated differential translation analysis
4. Compatibility with public databases storing processed Ribo-Seq data
