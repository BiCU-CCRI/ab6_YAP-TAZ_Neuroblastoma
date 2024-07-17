#Copywrite Valentina Boeva, 2017

# >>> SOURCE LICENSE >>>
#   This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation (www.fsf.org); either version 2 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# A copy of the GNU General Public License is available at
# http://www.fsf.org/licensing/licenses
# 
# >>> END OF LICENSE >>>

#################################################################################
# to run:
# cat runLILY.R | R --slave --args SAMPLE_NAME OUTPUT_DIR Distance_toStitch distFromTSS transcriptome_GFF_file faiFileToCreateBigWig
#################################################################################
Export_my_BED <- function(obj, path_to_export){
  BED <- data.frame(seqnames = seqnames(obj),
                    starts = start(obj)-1,
                    ends = end(obj),
                    strands = strand(obj))
  file_name <- deparse(substitute(obj))
  write.table(BED, file = file.path(path_to_export, paste0(file_name, ".bed")), quote=F, sep="\t", row.names=F, col.names=F)
}

dr_list <- list.dirs("/home/rstudio/workspace/neuroblastoma/CnR_LILY/20240604/", recursive = FALSE)

for(folder_name in dr_list){
  
  if(folder_name == "/home/rstudio/workspace/neuroblastoma/CnR_LILY/20240604//SK-N-SH-M_H3K27ac_R2"){next}
  
  args <- c(dummy1="",dummy2="",dummy3="",
            sampleName = paste0(folder_name,"/", basename(folder_name)),
            OUTPUT_DIR = paste0(folder_name, "/out"),
            maxDistanceToStitch = 12500,
            distFromTSS = 3000,
            transcriptomeFile = "/home/rstudio/workspace/neuroblastoma/resources/hg38_refseq_noChr.ucsc",
            faiFileToCreateBigWig = "/home/rstudio/workspace/neuroblastoma/resources/Homo_sapiens_Genome.GRCh38.102.fa.fai")
  
  sampleName=args[4]
  OUTPUT_DIR=args[5]
  maxDistanceToStitch=as.numeric(args[6])
  distFromTSS=as.numeric(args[7])
  transcriptomeFile=args[8]
  faiFileToCreateBigWig=args[9]
  
  cat (paste("..Working with sample",sampleName,"\n"))
  cat (paste("..Will stitch enhancers at maximal distance of",maxDistanceToStitch,"\n"))
  cat (paste("..Will consider promoter regions as regions around +-",distFromTSS,"from gene TSS\n"))
  
  
  #################check that all files and directories exist ####################
  
  dir.create(OUTPUT_DIR, showWarnings = FALSE)
  setwd(OUTPUT_DIR)
  cat (paste("..Will write the output into ",OUTPUT_DIR,"\n"))
  
  if(!file.exists(transcriptomeFile)) {
    cat ("Error:Please prodive a valid path to a file with transcriptome information\nFor example from here: https://github.com/linlabbcm/rose2/tree/master/rose2/annotation\nEXIT!\n")
  }
  
  if(!file.exists(paste0(sampleName,"_peaks.narrowPeak"))) {
    cat ("Error:Please prodive a valid path to output files of HMCan\n")
  }
  
  
  #####################################################
  
  suppressMessages(library(rtracklayer))
  
  ################# FUNCTIONS #########################
  
  
  import.bed3<- function(filename){
    peaks=read.table(filename)
    return(GRanges(peaks[,1],IRanges(peaks[,2],peaks[,3])))
  }
  
  import.bed6<- function(filename){
    peaks=read.table(filename)
    return(GRanges(peaks[,1],IRanges(peaks[,2],peaks[,3]),score=peaks[,5]))
  }
  
  import.ucsc<- function(filename){
    peaks=read.table(filename,header = T,comment.char = "!")
    peaks=GRanges(peaks$chrom,IRanges(peaks$txStart,peaks$txEnd),strand = peaks$strand)
    return(peaks)
  }
  
  
  getThresholdOnPeakScore <- function (filename) {
    dataTable <-read.table(filename, header=F);  
    lengths = dataTable$V3-dataTable$V2  
    lengths = (lengths - 1)/50 + 1  
    scores=dataTable$V5
    
    tranch=200
    ascendingScores = rev(scores)
    lengthForAscendingScores = rev(lengths)
    
    x = NULL
    y = NULL
    z=NULL
    for (i in c(0:(floor(length(lengths)/tranch)-1))) {
      tt=c((tranch*i+1):(tranch*i+tranch))
      m=mean(ascendingScores[tt]/lengthForAscendingScores[tt])
      s=sd(ascendingScores[tt]/lengthForAscendingScores[tt])
      x=c(x,max(ascendingScores[tt]))
      y=c(y,m)
      z=c(z,s)
    }
    
    thresholdToReturn=x[min(which(y>=min(y[which(x>5)])))]
    if (thresholdToReturn<2) {thresholdToReturn=2} #warning!!!
    return(thresholdToReturn)
  }
  
  output_peaks_density <- function(enhancersStitched,enhancers,promoters,bwFile,outputFile){
    
    bedRegions = enhancersStitched
    strand(bedRegions)="*"; 
    bedRegions$score=0; 
    densities = import.bw(bwFile)  
    gc()
    overlaps=findOverlaps(bedRegions,densities)  
    ensids <- densities$score[subjectHits(overlaps)]
    x <- tapply(ensids, queryHits(overlaps), sum)
    bedRegions$score[unique(queryHits(overlaps))]=x
    rm(overlaps)
    rm(densities)
    gc()  
    strand(bedRegions)="+"; 
    bedRegions=bedRegions[order(bedRegions$score,decreasing=T)]
    cutoff=calculate_cutoff(bedRegions$score,F)
    bedRegions$name="enhancer"
    bedRegions$name[which(bedRegions$score>cutoff)]="SE"
    
    #unstitch enhancers:
    strand(enhancers)="+";enhancers$score=0;enhancers$name="enhancer"; 
    tt1=which(countOverlaps(enhancers,bedRegions[which(bedRegions$name=="enhancer")])>0)
    tt2=which(countOverlaps(bedRegions,enhancers)>0 & bedRegions$name=="enhancer")
    
    bedRegions=bedRegions[-tt2]
    bedRegions=c(bedRegions,enhancers[tt1])
    
    #find promoters intersecting with enhancers (not SEs) and call these enhancers: promoters
    promoters$score=0
    promoters$name="promoter"
    strand(promoters)="+"
    
    enh_promoters=c(promoters,bedRegions[which(bedRegions$name=="enhancer")])
    enh_promoters=resize(enh_promoters,20+width(enh_promoters))
    enh_promoters=reduce(enh_promoters)
    enh_promoters=resize(enh_promoters,-20+width(enh_promoters))
    
    tt1=which(countOverlaps(enh_promoters,promoters)>0)
    tt2=which(countOverlaps(enh_promoters,promoters)==0)
    largePromoters=enh_promoters[tt1]
    largeEnh=enh_promoters[tt2]
    
    largePromoters$score=0;largePromoters$name="promoter"
    largeEnh$score=0;largeEnh$name="enhancer"
    tt1=which(countOverlaps(largeEnh,bedRegions[which(bedRegions$name=="enhancer")])>0)
    tt2=which(countOverlaps(bedRegions,largeEnh)>0 & bedRegions$name=="enhancer")
    bedRegions=bedRegions[-tt2]
    bedRegions=c(bedRegions,largeEnh[tt1])
    
    finalSet=c(bedRegions,largePromoters)
    strand(finalSet)="*"; 
    finalSet$score=0; 
    densities = import.bw(bwFile)  
    gc()
    overlaps=findOverlaps(finalSet,densities)  
    ensids <- densities$score[subjectHits(overlaps)]
    x <- tapply(ensids, queryHits(overlaps), sum)
    finalSet$score[unique(queryHits(overlaps))]=x
    rm(overlaps)
    rm(densities)
    gc()  
    strand(finalSet)="+"; 
    
    export.bed(finalSet,outputFile)
  }
  
  #============================================================================
  #==============SUPER-ENHANCER CALLING AND PLOTTING FUNCTIONS=================
  #======   from the original ROSE package developed by Charles Lin  (c) ======
  #======   ROSE licence:    http://younglab.wi.mit.edu/ROSE/LICENSE.txt ======
  #============================================================================
  
  #this is an accessory function, that determines the number of points below a diagnoal passing through [x,yPt]
  numPts_below_line <- function(myVector,slope,x){
    yPt <- myVector[x]
    b <- yPt-(slope*x)
    xPts <- 1:length(myVector)
    return(sum(myVector<=(xPts*slope+b)))
  }
  
  #This function calculates the cutoff by sliding a diagonal line and finding where it is tangential (or as close as possible)
  calculate_cutoff <- function(inputVector, drawPlot=FALSE,...){
    print("this version will try to get more than 600 SEs")
    t1=600
    
    inputVector <- sort(inputVector)
    inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
    
    numberOfSE=0
    while (numberOfSE<t1) {
      slope <- (max(inputVector)-min(inputVector))/(length(inputVector)) #This is the slope of the line we want to slide. This is the diagonal.
      
      xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
      y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
      
      numberOfSE = length(which(inputVector>y_cutoff))
      inputVector=inputVector[-length(inputVector)]
    }
    
    if(drawPlot){  #if TRUE, draw the plot
      plot(1:length(inputVector), inputVector,type="l",...)
      b <- y_cutoff-(slope* xPt)
      abline(v= xPt,h= y_cutoff,lty=2,col=8)
      points(xPt,y_cutoff,pch=16,cex=0.9,col=2)
      abline(coef=c(b,slope),col=2)
      title(paste("x=",xPt,"\ny=",signif(y_cutoff,3),"\nFold over Median=",signif(y_cutoff/median(inputVector),3),"x\nFold over Mean=",signif(y_cutoff/mean(inputVector),3),"x",sep=""))
      axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="pink",col="pink") #Number of regions with zero signal
    }
    return(y_cutoff)
  }
  
  ######## ######## ######## ######## ######## ######## ######## ######## ######## 
  ######## select threshold based on Narrow peaks:
  
  threshold_score = getThresholdOnPeakScore(paste0(sampleName,"_peaks.narrowPeak"))
  cat (paste("..Minimal peak score set to ",threshold_score,"\n"))
  
  ####### read regions and select regions with high score:
  
  regionsToStitch=import.bed6(paste0(sampleName,"_regions.bed"))
  tt=which(regionsToStitch$score>=threshold_score)
  cat (paste("will select ",length(tt), " regions out of ",length(regionsToStitch)," initial peaks\n"))
  
  regionsToStitch=regionsToStitch[tt]
  
  ######## read gene info and extract TSS:
  if (substr(transcriptomeFile,start = nchar(transcriptomeFile)-2,nchar(transcriptomeFile))=="gff") {
    transcripts=import(transcriptomeFile)
  }
  if  (substr(transcriptomeFile,start = nchar(transcriptomeFile)-3,nchar(transcriptomeFile))=="ucsc"){
    transcripts=import.ucsc(transcriptomeFile)
  }
  
  cat(paste("..Read file",transcriptomeFile,"\n"))
  
  TSSRegions=GRanges(chrom(transcripts),IRanges(start(transcripts)-distFromTSS,start(transcripts)+distFromTSS-1))
  tt=which(strand(transcripts)=='-')
  TSSRegions[tt]=GRanges(chrom(transcripts[tt]),IRanges(end(transcripts[tt])-distFromTSS,end(transcripts[tt])+distFromTSS-1))
  
  cat(paste("..Created",length(TSSRegions),"TSS regions\n"))
  
  ########## remove TSS peaks from enhancers:
  enhancers=setdiff(regionsToStitch,TSSRegions,ignore.strand=T)
  promoters=intersect(regionsToStitch,TSSRegions,ignore.strand=T)
  
  enhancers=IRanges::setdiff(regionsToStitch,TSSRegions,ignore.strand=T)
  promoters=intersect(regionsToStitch,TSSRegions,ignore.strand=T)
  ########## stitch enhancers:
  enhancersStitched=resize(enhancers,maxDistanceToStitch+width(enhancers))
  enhancersStitched=reduce(enhancersStitched)
  enhancersStitched=resize(enhancersStitched,-maxDistanceToStitch+width(enhancersStitched))
  
  cat(paste("..Created",length(enhancersStitched),"stitched regions\n"))
  #export(enhancersStitched,paste0(sampleName,".stitched_regions.bed"))
  
  ######### create big wig file out of wig if .bw does not exists ######
  hasBW=F
  bwFile=paste0(sampleName,".wig.bw")
  if(file.exists(bwFile)) {
    hasBW = TRUE
  }
  
  if (!hasBW) {
    bwFile=paste0(sampleName,".bw")
    if(file.exists(bwFile)) {
      hasBW = TRUE
    }
  }
  
  if (!hasBW) {
    cat ("Will create bw file using 'wigToBigWig'\nCheck what wigToBigWig is added to your PATH!\nEx: PATH=$PATH:/usr/local/bin/ucsc_tools/wigToBigWig\n")
    outSys=system(paste("wigToBigWig -clip", paste0(sampleName,".wig"), faiFileToCreateBigWig, paste0(sampleName,".bw")))
    bwFile=paste0(sampleName,".bw")
    if (outSys==0) {
      cat (paste(bwFile,"was created\n"))
    }
    if(file.exists(bwFile)) {
      hasBW = TRUE
    }else {
      cat(paste0("Could not create bigwig file out of ",paste0(sampleName,".wig"),"\nPlease create this file yourself and rerun\n"))
      cat (paste0("Ex: wigToBigWig -clip ", paste0(sampleName,".wig"), " YOUR_PATH_TO/Human/hg19/hg19.fa.fai ", paste0(sampleName,".bw\n")))
      cat ("wigToBigWig can be downloaded here: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/\n")
      quit()
    }
  }
  
  outputFile=paste0(sampleName,".scores.bed")
  outputFile=basename(outputFile)
  cat (paste("..Printing SEs, enhancers and promoters with their scores into",outputFile,"\n"))
  output_peaks_density(enhancersStitched,enhancers,promoters,bwFile,outputFile)
}
setwd("~/workspace/")

library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(GenomicRanges)
library(magrittr)

import::from(
  .from = "~/workspace/neuroblastoma/resources/utilityScripts.R",
  "pheatmap.type"
)

# Compose a list of enhancers
for(folder_name in dr_list){
  file.remove(paste0(folder_name, "/out/core"))
}
SE_list <- list()
for(folder_name in dr_list){
  if(dir.exists(paste0(folder_name, "/out/"))){
    SE_list[[basename(folder_name)]] <- read.delim(file = paste0(folder_name, "/out/", basename(folder_name), ".scores.bed"),
                                                   header = FALSE)
  }
}
SE_list <- lapply(SE_list, function(x){
  x <- x %>% filter(V4 == "SE")
  return(GRanges(seqnames = x$V1,
                 ranges = IRanges(start = x$V2,
                                  end = x$V3,
                                  names = paste0(x$V4, "_", x$V5))))
}
)


#CLB-Ma Enhancers
Grange_CLB_Ma_A_merged <- reduce(c(SE_list[[1]], SE_list[[2]]), drop.empty.ranges = F, min.gapwidth = 1)
Grange_CLB_Ma_M_merged <- reduce(c(SE_list[[3]], SE_list[[4]]), drop.empty.ranges = F, min.gapwidth = 1)

Grange_CLB_Ma_AM_reduced <- findOverlapPairs(Grange_CLB_Ma_M_merged, Grange_CLB_Ma_A_merged)
Grange_CLB_Ma_AM_reduced <- reduce(c(Grange_CLB_Ma_AM_reduced@first, c(Grange_CLB_Ma_AM_reduced@second)))

Grange_CLB_Ma_A <- findOverlaps(Grange_CLB_Ma_A_merged, Grange_CLB_Ma_AM_reduced)
Grange_CLB_Ma_A <- Grange_CLB_Ma_A_merged[-Grange_CLB_Ma_A@from,]

Grange_CLB_Ma_M <- findOverlaps(Grange_CLB_Ma_M_merged, Grange_CLB_Ma_AM_reduced)
Grange_CLB_Ma_M <- Grange_CLB_Ma_M_merged[-Grange_CLB_Ma_M@from,]

# export the coordinates 
Export_my_BED(Grange_CLB_Ma_A, "~/workspace/temp_results/BEDs/")
Export_my_BED(Grange_CLB_Ma_M, "~/workspace/temp_results/BEDs/")
Export_my_BED(Grange_CLB_Ma_AM_reduced, "~/workspace/temp_results/BEDs/")


# Normalization to the length of a SE 
tmp <- data.frame(
  Jun = c(countOverlaps(Grange_CLB_Ma_A,  DBobj_list_cons$Jun)/width(Grange_CLB_Ma_A)*10^3,
          countOverlaps(Grange_CLB_Ma_M,  DBobj_list_cons$Jun)/width(Grange_CLB_Ma_M)*10^3,
          countOverlaps(Grange_CLB_Ma_AM_reduced,  DBobj_list_cons$Jun)/width(Grange_CLB_Ma_AM_reduced)*10^3
  ),
  TAZ = c(countOverlaps(Grange_CLB_Ma_A,  DBobj_list_cons$TAZ)/width(Grange_CLB_Ma_A)*10^3,
          countOverlaps(Grange_CLB_Ma_M,  DBobj_list_cons$TAZ)/width(Grange_CLB_Ma_M)*10^3,
          countOverlaps(Grange_CLB_Ma_AM_reduced,  DBobj_list_cons$TAZ)/width(Grange_CLB_Ma_AM_reduced)*10^3
  ),
  YAP = c(countOverlaps(Grange_CLB_Ma_A,  DBobj_list_cons$YAP)/width(Grange_CLB_Ma_A)*10^3,
          countOverlaps(Grange_CLB_Ma_M,  DBobj_list_cons$YAP)/width(Grange_CLB_Ma_M)*10^3,
          countOverlaps(Grange_CLB_Ma_AM_reduced,  DBobj_list_cons$YAP)/width(Grange_CLB_Ma_AM_reduced)*10^3
  )
)

tmp$group <- c(rep("CLB_Ma_A", length(Grange_CLB_Ma_A)),
               rep("CLB_Ma_M", length(Grange_CLB_Ma_M)),
               rep("CLB_Ma_AM", length(Grange_CLB_Ma_AM_reduced)))
tmp$id <- c(paste0("CLB_Ma_A_", 1:length(Grange_CLB_Ma_A)),
            paste0("CLB_Ma_M_", 1:length(Grange_CLB_Ma_M)),
            paste0("CLB_Ma_AM_", 1:length(Grange_CLB_Ma_AM_reduced)))

rownames(tmp) <- tmp$id
annRow <-  data.frame(group = tmp$group)
row.names(annRow) <- tmp$id
pheatmap.type(tmp[ ,1:3], 
              annRow = annRow, 
              type = "group",
              doTranspose=FALSE, 
              conditions="Auto",
              cluster_cols = FALSE, 
              show_rownames = FALSE,
              color = colorRampPalette((brewer.pal(n = 7, name = "Reds")))(20))

# Now - SKN-N
Grange_SK_N_SH_A_merged <- reduce(c(SE_list[[7]], SE_list[[8]]), drop.empty.ranges = F, min.gapwidth = 1)
Grange_SK_N_SH_M_merged <- reduce(c(SE_list[[9]]), drop.empty.ranges = F, min.gapwidth = 1)

Grange_SK_N_SH_AM_reduced <- findOverlapPairs(Grange_SK_N_SH_M_merged, Grange_SK_N_SH_A_merged)
Grange_SK_N_SH_AM_reduced <- reduce(c(Grange_SK_N_SH_AM_reduced@first, c(Grange_SK_N_SH_AM_reduced@second)))

Grange_SK_N_SH_A <- findOverlaps(Grange_SK_N_SH_A_merged, Grange_SK_N_SH_AM_reduced)
Grange_SK_N_SH_A <- Grange_SK_N_SH_A_merged[-Grange_SK_N_SH_A@from,]

Grange_SK_N_SH_M <- findOverlaps(Grange_SK_N_SH_M_merged, Grange_SK_N_SH_AM_reduced)
Grange_SK_N_SH_M <- Grange_SK_N_SH_M_merged[-Grange_SK_N_SH_M@from,]

Export_my_BED(Grange_SK_N_SH_A, "~/workspace/temp_results/BEDs/")
Export_my_BED(Grange_SK_N_SH_M, "~/workspace/temp_results/BEDs/")
Export_my_BED(Grange_SK_N_SH_AM_reduced, "~/workspace/temp_results/BEDs/")

tmp <- data.frame(
  Jun = c(countOverlaps(Grange_SK_N_SH_A,  DBobj_list_cons$Jun)/width(Grange_SK_N_SH_A)*10^3,
          countOverlaps(Grange_SK_N_SH_M,  DBobj_list_cons$Jun)/width(Grange_SK_N_SH_M)*10^3, 
          countOverlaps(Grange_SK_N_SH_AM_reduced,  DBobj_list_cons$Jun)/width(Grange_SK_N_SH_AM_reduced)*10^3
  ),
  TAZ = c(countOverlaps(Grange_SK_N_SH_A,  DBobj_list_cons$TAZ)/width(Grange_SK_N_SH_A)*10^3,
          countOverlaps(Grange_SK_N_SH_M,  DBobj_list_cons$TAZ)/width(Grange_SK_N_SH_M)*10^3,
          countOverlaps(Grange_SK_N_SH_AM_reduced,  DBobj_list_cons$TAZ)/width(Grange_SK_N_SH_AM_reduced)*10^3
  ),
  YAP = c(countOverlaps(Grange_SK_N_SH_A,  DBobj_list_cons$YAP)/width(Grange_SK_N_SH_A)*10^3,
          countOverlaps(Grange_SK_N_SH_M,  DBobj_list_cons$YAP)/width(Grange_SK_N_SH_M)*10^3,
          countOverlaps(Grange_SK_N_SH_AM_reduced,  DBobj_list_cons$YAP)/width(Grange_SK_N_SH_AM_reduced)*10^3
  )
)
tmp$group <- c(rep("SK_N_SH_A", length(Grange_SK_N_SH_A)),
               rep("SK_N_SH_M", length(Grange_SK_N_SH_M)),
               rep("SK_N_SH_AM", length(Grange_SK_N_SH_AM_reduced)))
tmp$id <- c(paste0("SK_N_SH_A_", 1:length(Grange_SK_N_SH_A)),
            paste0("SK_N_SH_M_", 1:length(Grange_SK_N_SH_M)),
            paste0("SK_N_SH_AM_", 1:length(Grange_SK_N_SH_AM_reduced)))

rownames(tmp) <- tmp$id
annRow <-  data.frame(group = tmp$group)
row.names(annRow) <- tmp$id
pheatmap.type(tmp[ ,1:3], 
              annRow = annRow, 
              type = "group",
              doTranspose=FALSE, 
              conditions="Auto", 
              cluster_cols = FALSE, 
              show_rownames = FALSE,
              color = colorRampPalette((brewer.pal(n = 7, name = "Reds")))(20))


#list of peaks
library(ChIPpeakAnno)

ol <- ChIPpeakAnno::findOverlapsOfPeaks(Grange_SK_N_SH_A_merged, Grange_CLB_Ma_A_merged)
ChIPpeakAnno::makeVennDiagram(ol,
                              fill = c("#f5c9b1", "#f1b5ab"), # circle fill color
                              col = c("black", "black"), # circle border color
                              cat.col = c("#D55E00", "#0072B2")
)

ol <- ChIPpeakAnno::findOverlapsOfPeaks(Grange_SK_N_SH_M_merged, Grange_CLB_Ma_M_merged)
ChIPpeakAnno::makeVennDiagram(ol,
                              fill = c("#f5c9b1", "#f1b5ab"), # circle fill color
                              col = c("black", "black"), # circle border color
                              cat.col = c("#D55E00", "#0072B2")
)


## try to use only SEs that are present in both CLB and SKN
Grange_CLB_Ma_A_merged <- reduce(c(SE_list[[1]], SE_list[[2]]), drop.empty.ranges = F, min.gapwidth = 1)
Grange_SK_N_SH_A_merged <- reduce(c(SE_list[[7]], SE_list[[8]]), drop.empty.ranges = F, min.gapwidth = 1)

Grange_CLB_Ma_M_merged <- reduce(c(SE_list[[3]], SE_list[[4]]), drop.empty.ranges = F, min.gapwidth = 1)
Grange_SK_N_SH_M_merged <- reduce(c(SE_list[[9]]), drop.empty.ranges = F, min.gapwidth = 1)

CLB_SKN_A_merged <- findOverlapPairs(Grange_SK_N_SH_A_merged, Grange_CLB_Ma_A_merged)
CLB_SKN_A_merged <- reduce(c(CLB_SKN_A_merged@first, c(CLB_SKN_A_merged@second)))

CLB_SKN_M_merged <- findOverlapPairs(Grange_SK_N_SH_M_merged, Grange_CLB_Ma_M_merged)
CLB_SKN_M_merged <- reduce(c(CLB_SKN_M_merged@first, c(CLB_SKN_M_merged@second)))
  
CLB_SKN_AM <- findOverlapPairs(CLB_SKN_A_merged, CLB_SKN_M_merged)
CLB_SKN_AM <- reduce(c(CLB_SKN_AM@first, c(CLB_SKN_AM@second)))

CLB_SKN_A <- findOverlaps(CLB_SKN_A_merged, CLB_SKN_AM)
CLB_SKN_A <- CLB_SKN_A_merged[-CLB_SKN_A@from,]

CLB_SKN_M <- findOverlaps(CLB_SKN_M_merged, CLB_SKN_AM)
CLB_SKN_M <- CLB_SKN_M_merged[-CLB_SKN_M@from,]

Export_my_BED(CLB_SKN_A, "~/workspace/temp_results/BEDs/")
Export_my_BED(CLB_SKN_M, "~/workspace/temp_results/BEDs/")
Export_my_BED(CLB_SKN_AM, "~/workspace/temp_results/BEDs/")

tmp <- data.frame(
  Jun = c(countOverlaps(CLB_SKN_A,  DBobj_list_cons$Jun)/width(CLB_SKN_A)*10^3, 
          countOverlaps(CLB_SKN_M,  DBobj_list_cons$Jun)/width(CLB_SKN_M)*10^3,
          countOverlaps(CLB_SKN_AM, DBobj_list_cons$Jun)/width(CLB_SKN_AM)*10^3
  ),
  TAZ = c(countOverlaps(CLB_SKN_A,  DBobj_list_cons$TAZ)/width(CLB_SKN_A)*10^3, 
          countOverlaps(CLB_SKN_M,  DBobj_list_cons$TAZ)/width(CLB_SKN_M)*10^3,
          countOverlaps(CLB_SKN_AM, DBobj_list_cons$TAZ)/width(CLB_SKN_AM)*10^3
  ),
  YAP = c(countOverlaps(CLB_SKN_A,  DBobj_list_cons$YAP)/width(CLB_SKN_A)*10^3, 
          countOverlaps(CLB_SKN_M,  DBobj_list_cons$YAP)/width(CLB_SKN_M)*10^3,
          countOverlaps(CLB_SKN_AM, DBobj_list_cons$YAP)/width(CLB_SKN_AM)*10^3
  )
)
tmp$group <- c(rep("CLB_SKN_A", length(CLB_SKN_A)),
               rep("CLB_SKN_M", length(CLB_SKN_M)),
               rep("CLB_SKN_AM", length(CLB_SKN_AM)))
tmp$id <- c(paste0("CLB_SKN_A_", 1:length(CLB_SKN_A)),
            paste0("CLB_SKN_M_", 1:length(CLB_SKN_M)),
            paste0("CLB_SKN_AM_", 1:length(CLB_SKN_AM)))

rownames(tmp) <- tmp$id
annRow <-  data.frame(group = tmp$group)
row.names(annRow) <- tmp$id
pheatmap.type(tmp[ ,1:3], 
              annRow = annRow, 
              type = "group",
              doTranspose=FALSE, 
              conditions="Auto", 
              cluster_cols = FALSE, 
              show_rownames = FALSE,
              color = colorRampPalette((brewer.pal(n = 7, name = "Reds")))(20))

# alternative - only super unique peaks
ol <- ChIPpeakAnno::findOverlapsOfPeaks(Grange_CLB_Ma_A_merged, 
                                        Grange_SK_N_SH_A_merged, 
                                        Grange_CLB_Ma_M_merged, 
                                        Grange_SK_N_SH_M_merged)
ChIPpeakAnno::makeVennDiagram(ol,
                              fill = c("navyblue", "lightblue", "red", "firebrick"), # circle fill color
                              col = c("black", "black", "black", "black"), # circle border color
                              cat.col = c( "black", "black", "black", "black")
)
CLB_SKN_A_u <- ol[["peaklist"]][["Grange_CLB_Ma_A_merged///Grange_SK_N_SH_A_merged"]]
CLB_SKN_M_u <- ol[["peaklist"]][["Grange_CLB_Ma_M_merged///Grange_SK_N_SH_M_merged"]]
CLB_SKN_AM_u <- ol[["peaklist"]][["Grange_CLB_Ma_A_merged///Grange_SK_N_SH_A_merged///Grange_CLB_Ma_M_merged///Grange_SK_N_SH_M_merged"]]

tmp <- data.frame(
  Jun = c(countOverlaps( CLB_SKN_A_u,  DBobj_list_cons$Jun)/width(CLB_SKN_A_u)*10^3, 
          countOverlaps( CLB_SKN_M_u,  DBobj_list_cons$Jun)/width(CLB_SKN_M_u)*10^3,
          countOverlaps(CLB_SKN_AM_u, DBobj_list_cons$Jun)/width(CLB_SKN_AM_u)*10^3
  ),
  TAZ = c(countOverlaps( CLB_SKN_A_u,  DBobj_list_cons$TAZ)/width(CLB_SKN_A_u)*10^3, 
          countOverlaps( CLB_SKN_M_u,  DBobj_list_cons$TAZ)/width(CLB_SKN_M_u)*10^3,
          countOverlaps(CLB_SKN_AM_u, DBobj_list_cons$TAZ)/width(CLB_SKN_AM_u)*10^3
  ),
  YAP = c(countOverlaps( CLB_SKN_A_u,  DBobj_list_cons$YAP)/width(CLB_SKN_A_u)*10^3, 
          countOverlaps( CLB_SKN_M_u,  DBobj_list_cons$YAP)/width(CLB_SKN_M_u)*10^3,
          countOverlaps(CLB_SKN_AM_u, DBobj_list_cons$YAP)/width(CLB_SKN_AM_u)*10^3
  )
)
tmp$group <- c(rep("CLB_SKN_A", length(CLB_SKN_A_u)),
               rep("CLB_SKN_M", length(CLB_SKN_M_u)),
               rep("CLB_SKN_AM", length(CLB_SKN_AM_u)))
tmp$id <- c(paste0("CLB_SKN_A_", 1:length(CLB_SKN_A_u)),
            paste0("CLB_SKN_M_", 1:length(CLB_SKN_M_u)),
            paste0("CLB_SKN_AM_", 1:length(CLB_SKN_AM_u)))

library(ggpubr)
library(rstatix)

rownames(tmp) <- tmp$id
annRow <-  data.frame(group = tmp$group)
row.names(annRow) <- tmp$id
pheatmap.type(tmp[ ,1:3], 
              annRow = annRow, 
              type = "group",
              doTranspose=FALSE, 
              conditions="Auto", 
              cluster_cols = FALSE, 
              show_rownames = FALSE,
              color = colorRampPalette((brewer.pal(n = 7, name = "Reds")))(20))

tmp_long <- tmp %>% 
  dplyr::select(-id) %>%
  tidyr::pivot_longer(cols = c("Jun", "TAZ", "YAP"), names_to = "TF", values_to = "enrich_score")
  

my_comparisons <- list( 
  c("CLB_SKN_A", "CLB_SKN_AM"), 
  c("CLB_SKN_A", "CLB_SKN_M"),
  c("CLB_SKN_AM", "CLB_SKN_M") 
)

stat.test <- tmp_long %>%
  group_by(TF) %>%
  wilcox_test(enrich_score ~ group) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test <- stat.test %>% add_xy_position(x = "group")

ggboxplot(tmp_long ,
          x = "group",
          y = "enrich_score",
          fill = "group",
          facet.by = "TF",
          order = c("CLB_SKN_A", "CLB_SKN_AM", "CLB_SKN_M")
) + 
  stat_pvalue_manual(stat.test, label = "p.adj") +
  rotate_x_text(angle = 45)

ggboxplot(tmp_long ,
          x = "group",
          y = "enrich_score",
          fill = "group",
          facet.by = "TF",
          order = c("CLB_SKN_A", "CLB_SKN_AM", "CLB_SKN_M")
) + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif") +
  rotate_x_text(angle = 45)


+

+ # Add pairwise comparisons p-value
  stat_compare_means()


ggplot(tmp_long) +
  aes(x = TF, y = enrich_score, fill = group) +
  geom_boxplot() +
  scale_fill_manual(
    values = c(CLB_SKN_A = "#04d45c",
               CLB_SKN_AM = "#fc948c",
               CLB_SKN_M = "#84b4fc")
  ) +
  theme_minimal() +
  facet_wrap(vars(group))



# Calculate the enrichment of SEs in MES/ADR genes identified from RNA-seq
library(chipenrich)
library(biomaRt)


# BEDs must be annotated first using homer:
# "~/workspace/neuroblastoma/05_SEs_annotation.sh"
# Then load the files:
CLB_SKN_A_SEs_annot <- read.csv("~/workspace/temp_results/BEDs/CLB_SKN_A_homer_annot.csv", sep = "\t")
CLB_SKN_M_SEs_annot <- read.csv("~/workspace/temp_results/BEDs/CLB_SKN_M_homer_annot.csv", sep = "\t")
CLB_SKN_AM_SEs_annot <- read.csv("~/workspace/temp_results/BEDs/CLB_SKN_AM_homer_annot.csv", sep = "\t")

# load genes from RNA-seq MES/ADR-specific genes
RNA_SEQ_data <- openxlsx2::read_xlsx("~/workspace/neuroblastoma/results/20240515/cell_type_MES_vs_ADR.xlsx", sheet = 1)
mes_adrn_gene_list <- RNA_SEQ_data %>%
  mutate(Term = if_else(log2FoldChange < 0, "ADRN", "MES"))
mes_adrn_gene_list <- list(
  MES = mes_adrn_gene_list[mes_adrn_gene_list$Term == "MES", "gene_symbol"],
  ADR = mes_adrn_gene_list[mes_adrn_gene_list$Term == "ADRN", "gene_symbol"]
)

# GSEA
CLB_SKN_A_SEs_annot <- hypeR::hypeR(signature = CLB_SKN_A_SEs_annot$Gene.Name, 
                        genesets = mes_adrn_gene_list, 
                        test="hypergeometric"
                        )
CLB_SKN_M_SEs_annot <- hypeR::hypeR(signature = CLB_SKN_M_SEs_annot$Gene.Name, 
                    genesets = mes_adrn_gene_list, 
                    test="hypergeometric"
)
CLB_SKN_AM_SEs_annot <- hypeR::hypeR(signature = CLB_SKN_AM_SEs_annot$Gene.Name, 
                    genesets = mes_adrn_gene_list, 
                    test="hypergeometric"
)

SEs_GSEA_df <- NULL
SEs_GSEA_df <- dplyr::bind_rows(
  c(sample = "CLB_SKN_A_SEs", 
    type = "ADR",
    fdr = as.numeric(CLB_SKN_A_SEs_annot$data["ADR", "fdr"]), 
    overlap = as.numeric(CLB_SKN_A_SEs_annot$data["ADR", "overlap"])),
  c(sample = "CLB_SKN_A_SEs",
    type = "MES",
    fdr = as.numeric(CLB_SKN_A_SEs_annot$data["MES", "fdr"]), 
    overlap = as.numeric(CLB_SKN_A_SEs_annot$data["MES", "overlap"])),
  c(sample = "CLB_SKN_AM_SEs", 
    type = "ADR",
    fdr = as.numeric(CLB_SKN_AM_SEs_annot$data["ADR", "fdr"]), 
    overlap = as.numeric(CLB_SKN_AM_SEs_annot$data["ADR", "overlap"])), 
  c(sample = "CLB_SKN_AM_SEs", 
    type = "MES",
    fdr = as.numeric(CLB_SKN_AM_SEs_annot$data["MES", "fdr"]), 
    overlap = as.numeric(CLB_SKN_AM_SEs_annot$data["MES", "overlap"])),
  c(sample = "CLB_SKN_M_SEs", 
    type = "ADR",
    fdr = as.numeric(CLB_SKN_M_SEs_annot$data["ADR", "fdr"]), 
    overlap = as.numeric(CLB_SKN_M_SEs_annot$data["ADR", "overlap"])), 
  c(sample = "CLB_SKN_M_SEs", 
    type = "MES",
    fdr = as.numeric(CLB_SKN_M_SEs_annot$data["MES", "fdr"]), 
    overlap = as.numeric(CLB_SKN_M_SEs_annot$data["MES", "overlap"]))
)

SEs_GSEA_df$fdr %<>% as.numeric
SEs_GSEA_df$overlap %<>% as.numeric
SEs_GSEA_df$fdr <- -log10(SEs_GSEA_df$fdr + 0.0000000001)

SEs_GSEA_df %>%
  ggplot(aes(x=sample, y = type, color = fdr, size = overlap)) + 
  geom_point() +
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') +
  theme(axis.ticks = element_blank()) +
  labs(color = "-log10(FDR)") +
  #scale_color_distiller(palette = "Reds", direction = 1)
  scale_color_gradient(low = "grey", high = "red")

