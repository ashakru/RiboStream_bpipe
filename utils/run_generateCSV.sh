#!/bin/bash

#source /etc/profile.d/modules.sh

#module purge !!!!!!!!!!!!!
#module load use.own

#module load python/2.7.5

home_dir=$PWD
prog_dir=$RIBOSTREAM_HOME

python $prog_dir/utils/generateCSV.py "$@"
