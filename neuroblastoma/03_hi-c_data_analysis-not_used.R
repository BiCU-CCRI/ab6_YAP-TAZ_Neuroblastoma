#take a reduced dataset from 
#http://kobic.kr/3div/hic
# MSC induced. 
# only bids around chr11:102112447 - YAP

library(dplyr)
library(tidyr)
library(IRanges)
import::from(IRanges, IRanges)
import::from(GenomicRanges,GRanges)

#A549 dataset
A549_hic_data_path <- file.path("~/workspace/neuroblastoma/data_public/Hi-C/A549.chr11.covnorm.gz")
A549_hic_data <- read.table(A549_hic_data_path, header = T)
A549_hic_data <- A549_hic_data %>% separate(frag1, c('CHR1', "start1", "end1"))  %>% separate(frag2, c('CHR2', "start2", "end2"))
A549_hic_data <- A549_hic_data %>% select(-one_of('CHR1', "CHR2"))

#102120000 - is the closest start1
keep <- A549_hic_data$start1 == A549_hic_data[which.min(abs(102112447-as.numeric(A549_hic_data$start1))),]$start1
A549_hic_data_YAP <- A549_hic_data[keep,]
A549_hic_data_YAP <- A549_hic_data_YAP %>% mutate(biased_removed = capture_res/exp_value_capture)
#i will use value of 2 - totally arbitrary.
A549_hic_data_YAP <- A549_hic_data_YAP[A549_hic_data_YAP$biased_removed > 2,]
write.csv2(A549_hic_data_YAP, file = "~/workspace/neuroblastoma/data_public/Hi-C/processed/A549_hic_data_YAP_1_capture_res.csv")


#SKN dataset
SKN_hic_data_path <- file.path("~/workspace/neuroblastoma/data_public/Hi-C/Sk-N-SH.chr11.covnorm.gz")
SKN_hic_data <- read.table(SKN_hic_data_path, header = T)
SKN_hic_data <- SKN_hic_data %>% separate(frag1, c('CHR1', "start1", "end1"))  %>% separate(frag2, c('CHR2', "start2", "end2"))
SKN_hic_data <- SKN_hic_data %>% select(-one_of('CHR1', "CHR2"))

#102120000 - is the closest start1
keep <- SKN_hic_data$start1 == SKN_hic_data[which.min(abs(102112447-as.numeric(SKN_hic_data$start1))),]$start1
SKN_hic_data_YAP <- SKN_hic_data[keep,]
SKN_hic_data_YAP <- SKN_hic_data_YAP %>% mutate(biased_removed = capture_res/exp_value_capture)
#i will use value of 2 - totally arbitrary.
SKN_hic_data_YAP <- SKN_hic_data_YAP[SKN_hic_data_YAP$biased_removed > 2,]
write.csv2(SKN_hic_data_YAP, file = "~/workspace/neuroblastoma/data_public/Hi-C/processed/SKN_hic_data_YAP_1_capture_res.csv")


#SW dataset
SW_hic_data_path <- file.path("~/workspace/neuroblastoma/data_public/Hi-C/SW480_neg_control.chr11.covnorm.gz")
SW_hic_data <- read.table(SW_hic_data_path, header = T)
SW_hic_data <- SW_hic_data %>% separate(frag1, c('CHR1', "start1", "end1"))  %>% separate(frag2, c('CHR2', "start2", "end2"))
SW_hic_data <- SW_hic_data %>% select(-one_of('CHR1', "CHR2"))

#102120000 - is the closest start1
keep <- SW_hic_data$start1 == SW_hic_data[which.min(abs(102112447-as.numeric(SW_hic_data$start1))),]$start1
SW_hic_data_YAP <- SW_hic_data[keep,]
SW_hic_data_YAP <- SW_hic_data_YAP %>% mutate(biased_removed = capture_res/exp_value_capture)
#i will use value of 2 - totally arbitrary.
SW_hic_data_YAP <- SW_hic_data_YAP[SW_hic_data_YAP$biased_removed > 2,]
write.csv2(SW_hic_data_YAP, file = "~/workspace/neuroblastoma/data_public/Hi-C/processed/SW_hic_data_YAP_1_capture_res.csv")


#U2 dataset
U2_hic_data_path <- file.path("~/workspace/neuroblastoma/data_public/Hi-C/U2OS_G1.chr11.covnorm.gz")
U2_hic_data <- read.table(U2_hic_data_path, header = T)
U2_hic_data <- U2_hic_data %>% separate(frag1, c('CHR1', "start1", "end1"))  %>% separate(frag2, c('CHR2', "start2", "end2"))
U2_hic_data <- U2_hic_data %>% select(-one_of('CHR1', "CHR2"))

#102120000 - is the closest start1
keep <- U2_hic_data$start1 == U2_hic_data[which.min(abs(102112447-as.numeric(U2_hic_data$start1))),]$start1
U2_hic_data_YAP <- U2_hic_data[keep,]
U2_hic_data_YAP <- U2_hic_data_YAP %>% mutate(biased_removed = capture_res/exp_value_capture)
#i will use value of 2 - totally arbitrary.
U2_hic_data_YAP <- U2_hic_data_YAP[U2_hic_data_YAP$biased_removed > 2,]
write.csv2(U2_hic_data_YAP, file = "~/workspace/neuroblastoma/data_public/Hi-C/processed/U2_hic_data_YAP_1_capture_res.csv")


#WI dataset
WI_hic_data_path <- file.path("~/workspace/neuroblastoma/data_public/Hi-C/WI38_proliferative.chr11.covnorm.gz")
WI_hic_data <- read.table(WI_hic_data_path, header = T)
WI_hic_data <- WI_hic_data %>% separate(frag1, c('CHR1', "start1", "end1"))  %>% separate(frag2, c('CHR2', "start2", "end2"))
WI_hic_data <- WI_hic_data %>% select(-one_of('CHR1', "CHR2"))

#102120000 - is the closest start1
keep <- WI_hic_data$start1 == WI_hic_data[which.min(abs(102112447-as.numeric(WI_hic_data$start1))),]$start1
WI_hic_data_YAP <- WI_hic_data[keep,]
WI_hic_data_YAP <- WI_hic_data_YAP %>% mutate(biased_removed = capture_res/exp_value_capture)
#i will use value of 2 - totally arbitrary.
WI_hic_data_YAP <- WI_hic_data_YAP[WI_hic_data_YAP$biased_removed > 2,]
write.csv2(WI_hic_data_YAP, file = "~/workspace/neuroblastoma/data_public/Hi-C/processed/WI_hic_data_YAP_1_capture_res.csv")


#heladataset
HeLa_hic_data_path <- file.path("~/workspace/neuroblastoma/data_public/Hi-C/HeLaKyoto.chr11.covnorm.gz")
HeLa_hic_data <- read.table(HeLa_hic_data_path, header = T)
HeLa_hic_data <- HeLa_hic_data %>% separate(frag1, c('CHR1', "start1", "end1"))  %>% separate(frag2, c('CHR2', "start2", "end2"))
HeLa_hic_data <- HeLa_hic_data %>% select(-one_of('CHR1', "CHR2"))

#102120000 - is the closest start1
keep <- HeLa_hic_data$start1 == HeLa_hic_data[which.min(abs(102112447-as.numeric(HeLa_hic_data$start1))),]$start1
HeLa_hic_data_YAP <- HeLa_hic_data[keep,]
HeLa_hic_data_YAP <- HeLa_hic_data_YAP %>% mutate(biased_removed = capture_res/exp_value_capture)
#i will use value of 2 - totally arbitrary.
HeLa_hic_data_YAP <- HeLa_hic_data_YAP[HeLa_hic_data_YAP$biased_removed > 2,]
write.csv2(HeLa_hic_data_YAP, file = "~/workspace/neuroblastoma/data_public/Hi-C/processed/HeLa_hic_data_YAP_1_capture_res.csv")


#Now - scale2 format experiments

#fibroC dataset
fibroC_hic_data_path <- file.path("~/workspace/neuroblastoma/data_public/Hi-C/fibro_control.chr11.distnorm.scale2.gz")
fibroC_hic_data <- read.table(fibroC_hic_data_path, header = T) 
fibroC_ranges <- IRanges(names = fibroC_hic_data$idx, start = fibroC_hic_data$bin1, end = fibroC_hic_data$bin1+5000)

YAP_TSS <- IRanges(start = 102112447-3000, end = 102112447+1000)
keep <- subsetByOverlaps(fibroC_ranges, YAP_TSS)
fibroC_hic_data_YAP <- fibroC_hic_data[fibroC_hic_data$bin1 %in% start(keep),]
#i will use value of 2 - totally arbitrary.
fibroC_hic_data_YAP <- fibroC_hic_data_YAP[fibroC_hic_data_YAP$dist_foldchange > 2, ]
write.csv2(fibroC_hic_data_YAP, file = "~/workspace/neuroblastoma/data_public/Hi-C/processed/fibroC_hic_data_YAP_1_capture_res.csv")


#IMR90 dataset
IMR90_hic_data_path <- file.path("~/workspace/neuroblastoma/data_public/Hi-C/IMR90_fibroblast.chr11.distnorm.scale2.gz")
IMR90_hic_data <- read.table(IMR90_hic_data_path, header = T)
IMR90_ranges <- IRanges(names = IMR90_hic_data$idx, start = IMR90_hic_data$bin1, end = IMR90_hic_data$bin1+5000)

YAP_TSS <- IRanges(start = 102112447-3000, end = 102112447+1000)
keep <- subsetByOverlaps(IMR90_ranges, YAP_TSS)
IMR90_hic_data_YAP <- IMR90_hic_data[IMR90_hic_data$bin1 %in% start(keep),]
#i will use value of 2 - totally arbitrary.
IMR90_hic_data_YAP <- IMR90_hic_data_YAP[IMR90_hic_data_YAP$dist_foldchange > 2, ]
write.csv2(IMR90_hic_data_YAP, file = "~/workspace/neuroblastoma/data_public/Hi-C/processed/IMR90_hic_data_YAP_1_capture_res.csv")

#MCF dataset # THE DATASET IS TOO SMALL
MCF_hic_data_path <- file.path("~/workspace/neuroblastoma/data_public/Hi-C/MCF-10A_control.chr11.distnorm.scale2.gz")
MCF_hic_data <- read.table(MCF_hic_data_path, header = T)
keep <- MCF_hic_data$bin1 == MCF_hic_data[which.min(abs(102112447-as.numeric(MCF_hic_data$bin1))),]$bin1
MCF_hic_data_YAP <- MCF_hic_data[keep,]
#i will use value of 2 - totally arbitrary.
MCF_hic_data_YAP <- MCF_hic_data_YAP[MCF_hic_data_YAP$dist_foldchange > 2, ]
write.csv2(MCF_hic_data_YAP, file = "~/workspace/neuroblastoma/data_public/Hi-C/processed/MCF_hic_data_YAP_1_capture_res.csv")


#MSC dataset
MSC_hic_data_path <- file.path("~/workspace/neuroblastoma/data_public/Hi-C/MSC.chr11.distnorm.scale2.gz")
MSC_hic_data <- read.table(MSC_hic_data_path, header = T)
MSC_ranges <- IRanges(names = MSC_hic_data$idx, start = MSC_hic_data$bin1, end = MSC_hic_data$bin1+5000)
YAP_TSS <- IRanges(start = 102112447-3000, end = 102112447+1000)
keep <- subsetByOverlaps(MSC_ranges, YAP_TSS)

MSC_hic_data_YAP <- MSC_hic_data[MSC_hic_data$bin1 %in% start(keep),]
#i will use value of 2 - totally arbitrary.
MSC_hic_data_YAP <- MSC_hic_data_YAP[MSC_hic_data_YAP$dist_foldchange > 2, ]
write.csv2(MSC_hic_data_YAP, file = "~/workspace/neuroblastoma/data_public/Hi-C/processed/MSC_hic_data_YAP_1_capture_res.csv")



# extract scale table bins that are present at least 2 times
scale2_tbl <- table(rbind(IMR90_hic_data_YAP,MSC_hic_data_YAP, fibroC_hic_data_YAP)$bin2)
names(scale2_tbl[scale2_tbl>1])

scale2_bed <- data.frame(
chr = 'chr11',
start = as.numeric(names(scale2_tbl[scale2_tbl>1])),
end = as.numeric(names(scale2_tbl[scale2_tbl>1]))+5000)

#export as BED file 
write.csv2(scale2_bed, file = "~/workspace/neuroblastoma/data_public/Hi-C/processed/scale2_bed.csv")






##############
#
#Processing the full HiC dataset - not real need to do that.
#

file_list <- list.files(path = "~/workspace/neuroblastoma/data_public/Hi-C/full_data/MSC/", pattern = "*.gz", all.files = FALSE,
           full.names = T, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

for (file in file_list){
  tbl_file <- read.table(file)
  
}

##############
#
#Processing data from: https://sci-hub.ru/https://www.nature.com/articles/ncb3216
#

MDA_enhancers <- read.table("~/workspace/neuroblastoma/resources/MDA-MB-231_enhancers.csv", header=TRUE, sep = ";")
MDA_enhancers_IRanges <- GRanges(MDA_enhancers)

ATAC_dds <- readRDS("~/workspace/neuroblastoma/RDSs/ATAC_dds.RDS")

intersect(rowRanges(ATAC_dds), MDA_enhancers_IRanges)


