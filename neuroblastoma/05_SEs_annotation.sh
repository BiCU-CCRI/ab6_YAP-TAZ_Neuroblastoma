#!/bin/bash

annotatePeaks.pl ../temp_results/BEDs/CLB_SKN_A_chr.bed hg38 >../temp_results/BEDs/CLB_SKN_A_homer_annot.csv
annotatePeaks.pl ../temp_results/BEDs/CLB_SKN_M_chr.bed hg38 >../temp_results/BEDs/CLB_SKN_M_homer_annot.csv
annotatePeaks.pl ../temp_results/BEDs/CLB_SKN_AM_chr.bed hg38 >../temp_results/BEDs/CLB_SKN_AM_homer_annot.csv
