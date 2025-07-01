#!/bin/bash

# install HOMER
# http://homer.ucsd.edu/homer/

homer_dir=/home/rstudio/workspace/neuroblastoma/homer

mkdir ${homer_dir}
wget http://homer.ucsd.edu/homer/configureHomer.pl -O ${homer_dir}/configureHomer.pl
perl ${homer_dir}/configureHomer.pl -install
perl ${homer_dir}/configureHomer.pl -install hg38

PATH=$PATH:/home/rstudio/workspace/neuroblastoma/homer/bin/
