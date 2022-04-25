#!/bin/bash

#Script to perform initial preprocessing of provided fastq files

source /etc/profile.d/modules.sh

home_dir=$PWD
prog_dir=$RIBOSTREAM_HOME

echo $prog_dir
echo $home_dir

#Check if a directory with fastq files exists
dir="fastq"
if [[ ! -e $dir ]]; then
    mkdir $dir
elif [[ ! -d $dir ]]; then
    echo "$dir already exists but is not a directory" 1>&2
fi

#Check if configs directory (with all setting for each tool used)
dir="configs"
if [[ ! -e $dir ]]; then
    mkdir $dir
elif [[ ! -d $dir ]]; then
    echo "$dir already exists but is not a directory" 1>&2
fi

# Script that generates JSON files with setting for each sample
Rscript $prog_dir/utils/generateJSON.R $1 "_preproc"

# Save custer settings and tools parameters
cp $prog_dir/bpipe.config configs/
cp $home_dir/config.ini configs/

# Run the pipeline
bpipe run -n 5 -r $prog_dir/init_preproc.pipeline configs/*.json
