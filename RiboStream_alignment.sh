#!/bin/bash

# $1: file name of sample sheet (.csv)

source /etc/profile.d/modules.sh
home_dir=$PWD
prog_dir=$RIBOSTREAM_HOME

# TODO: for loop for folders

dir="fastqc"
if [[ ! -e $dir ]]; then
    mkdir $dir
elif [[ ! -d $dir ]]; then
    echo "$dir already exists but is not a directory" 1>&2
fi

dir="bam"
if [[ ! -e $dir ]]; then
    mkdir $dir
elif [[ ! -d $dir ]]; then
    echo "$dir already exists but is not a directory" 1>&2
fi

dir="RNA_transcriptome"
if [[ ! -e $dir ]]; then
    mkdir $dir
elif [[ ! -d $dir ]]; then
    echo "$dir already exists but is not a directory" 1>&2
fi

dir="configs"
if [[ ! -e $dir ]]; then
    mkdir $dir
elif [[ ! -d $dir ]]; then
    echo "$dir already exists but is not a directory" 1>&2
fi

# #generate sample sheet for each step
# module load gcc/4.9.2
# module load R/3.2.2

Rscript $prog_dir/utils/generateJSON.R $1 "_alignment"
Rscript $prog_dir/utils/generateExpJSON.R $1 "_experiment"
#Rscript $prog_dir/generateCondJSON.R $1 "_condition"
#sample_sheet=$1

# rename sampleSheet_compare:
# mv -u $2 sampleSheet_compare.csv

# save version the number of each tool used in pipeline
cp $prog_dir/bpipe.config configs/

#run the pipeline
bpipe run -r $prog_dir/init_alignment.pipeline configs/*
