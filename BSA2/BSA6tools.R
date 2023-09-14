#This is a function repository for the BSA experiment.
# Access by source("/Users/Daniel/Google_Drive/WashU/DNA/tools/BSAtools.R")

library(knitr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(reshape2)
library(tidyverse)
#library(QTLseqr)

##########################################################################################################################################
#' This function extracts info from an imported vcf for 1 sample, outputting a nice table.
#' 
#' @param VCF 
#' @param i integer. Samples in a vcf file start in the 10th column. Tabler takes the information from the (9+i)th column from the vcf.
#' @param minDepth Numeric. The minimal depth required at each position to call a different allele. Default=5.
#' @param cutoff Numeric. Minimum ratio of reads (as a number between 0 and 1) that supports either allele required to call an allele. Default=0.75. 
#' @param filtering Logical. If True, will fill up the column "filtGenotype". This is specially useful for single strains, not so much for BSA. The function is much faster if turned off. Default=F. 
#' @return a table

tabler=function(VCF, i=1, minDepth=5, cutoff=0.75, filtering=F, cutoff50=0.1){
  c=(stringr::str_split_fixed(VCF[,9+i],":",n=8)) #This splits the FORMAT itens from the vcf into columns.
  colnames(c)=c("Genotype_GT","Depth_DP","AD","Reference_RO","QR","Alternative_AO","QA","GL")
  #The two lines bellow convert any "." in the object c (columns 3, 8 and the rest, respectively) to "0,0".
  c[,3]=sapply(c[,3], USE.NAMES = F, FUN = function(x){if(x=="."){x="0,0"} else {x=x}  })
  c[,8]=sapply(c[,8], USE.NAMES = F, FUN = function(x){if(x=="."){x="0,0"} else {x=x}  })
  c=apply(c,c(1,2), FUN = function(x){if(x=="."){x=0} else {x=x}  }  )
  
  #This gets the the number of Reference reads and ALternative reads from the column AD from the object c.
  d=stringr::str_split_fixed(c[,3], ",",n=3) #This splits the two alternative allele counts into different columns.
  
  #Convert c and d to numeric
  c=apply(c,c(1,2),as.numeric) #Columns AD and GL will be converted to NA, generating errors. But that's ok, we've got what we needed from them. 
  d=apply(d,c(1,2),as.numeric)
  
  #This gets the possible alternative alleles. The split only happens when we have a second alternative allele.
  e=stringr::str_split_fixed(VCF[,5], ",",n=3)

  data2=data.frame(CHR=VCF[,1], POS=VCF[,2], Ref_allele= VCF[,4], ALT1=e[,1],ALT2=e[,2],ALT3=e[,3], Genotype=as.numeric(c[,1]),
                   Depth=as.numeric(c[,2]),Reference= as.numeric(c[,4]), Alternative1=d[,2], Alternative2=d[,3], 
                   filtGenotype=rep(".",nrow(e)), stringsAsFactors = FALSE) 
  data2$Alternative2=sapply(data2$Alternative2, USE.NAMES = F, FUN = function(x){if(is.na(x)){x=as.numeric(0)} else {x=x}  })
  
  data2$filtGenotype=as.character(data2$filtGenotype)
  data2=dplyr::mutate(data2, realDepth=Reference+Alternative1+Alternative2,
         Ref_percentage=Reference/realDepth,
         Alt1_percentage=Alternative1/realDepth,
         Alt2_percentage=Alternative2/realDepth)
  if(filtering==T){
  for(j in 1:nrow(data2)){ ## Loop for defining the filtered Genotype. Filter by minimum depth and Genotype cutoff. The cutoff is the
    # minimum percentage of reads necessary to define a genotype. Default = 0.75.
    #the cutoff50 is the distance to 0.5 to a position so it can be called as "Mixed", as in mixed strains.
    data2$filtGenotype
    if(is.na(data2$Reference[j])){x="noReads"}
    else if(data2$realDepth[j]< (minDepth+1)){data2[j,12]="lowDepth"}
    else if(data2$Reference[j]>cutoff*data2$realDepth[j]){data2[j,12]="Reference"}
    else if(data2$Alternative1[j]>cutoff*data2$realDepth[j]){data2[j,12]="Alternative1"}
    else if(data2$Alternative2[j]>cutoff*data2$realDepth[j]){data2[j,12]="Alternative2"}
    else if((data2$Ref_percentage[j]>(0.5-cutoff50) &  data2$Ref_percentage[j]<(0.5+cutoff50) ) &
            (data2$Alt1_percentage[j]>(0.5-cutoff50) &  data2$Alt1_percentage[j]<(0.5+cutoff50))  ){data2[j,12]="Mixed"}
    else if((data2$Alternative2[j]<=cutoff*data2$realDepth[j] &  data2$Alternative1[j]<=cutoff*data2$realDepth[j] &
            data2$Reference[j]<=cutoff*data2$realDepth[j])  ){data2[j,12]="undefinedGenotype"}
    
  }
}
  rm(c)
  rm(d)
  return(data2)
}

##########################

#' This function creates a list with each element being the comparison between two bulks. This will be in the format used by QTLseq package.
#' 
#' @param SNPset  List. The main input for this function is a list in which all   
#' @param lowBulk Character. The lowbulk of the comparison made.
#' @param highBulk Character. The highbulk of the comparison made.
#' @return data-frame that can be used as an element of a list.


picker=function(SNPset, lowBulk, highBulk){
  if(!lowBulk %in% names(SNPset) | !highBulk %in% names(SNPset)){
    stop("Confirm if the bulks' names are present in the SNP set (vcf file) provided.")
  }
  metaTable=data.frame(CHROM=SNPset[[lowBulk]]$CHR, POS=SNPset[[lowBulk]]$POS, REF=SNPset[[lowBulk]]$Ref_allele, 
                       ALT=SNPset[[lowBulk]]$ALT1, 
                       AD_REF.LOW=SNPset[[lowBulk]]$Reference,
                       AD_ALT.LOW=(SNPset[[lowBulk]]$Alternative1+(SNPset[[lowBulk]]$Alternative2)),
                       DP.LOW=SNPset[[lowBulk]]$realDepth,
                       GQ.LOW=rep(NA, nrow(SNPset[[lowBulk]])) ,
                       PL.LOW=rep(NA, nrow(SNPset[[lowBulk]])),
                       SNPindex.LOW=(SNPset[[lowBulk]]$Alternative1+(SNPset[[lowBulk]]$Alternative2))/SNPset[[lowBulk]]$realDepth,
                       AD_REF.HIGH=SNPset[[highBulk]]$Reference,
                       AD_ALT.HIGH=(SNPset[[highBulk]]$Alternative1+(SNPset[[highBulk]]$Alternative2)),
                       DP.HIGH=SNPset[[highBulk]]$realDepth,
                       GQ.HIGH=rep(NA, nrow(SNPset[[highBulk]])) ,
                       PL.HIGH=rep(NA, nrow(SNPset[[highBulk]])),
                       SNPindex.HIGH=(SNPset[[highBulk]]$Alternative1+(SNPset[[highBulk]]$Alternative2))/SNPset[[highBulk]]$realDepth, #AD_ALT.HIGH/DP.HIGH,
                       REF_FRQ=((SNPset[[highBulk]]$Reference+SNPset[[lowBulk]]$Reference)/(SNPset[[highBulk]]$realDepth+SNPset[[lowBulk]]$realDepth))
  )
  metaTable=dplyr::mutate(metaTable, deltaSNP=SNPindex.HIGH-SNPindex.LOW)
  return(metaTable)
}

##########################

#' This function creates a list with each element being the comparison between two bulks. This will be in the format used by QTLseq package.
#' Picker2 has some adaptations from picker 1 for the tabler function in awk. RealDepth instead of realDepth, and no Alt2 allele.
#' @param SNPset  List. The main input for this function is a list in which all   
#' @param lowBulk Character. The lowbulk of the comparison made.
#' @param highBulk Character. The highbulk of the comparison made.
#' @return data-frame that can be used as an element of a list.


picker2=function(SNPset, lowBulk, highBulk){
  if(!lowBulk %in% names(SNPset) | !highBulk %in% names(SNPset)){
    stop("Confirm if the bulks' names are present in the SNP set (vcf file) provided.")
  }
  
  metaTable=data.frame(CHROM=SNPset[[lowBulk]]$CHR, 
                       POS=SNPset[[lowBulk]]$POS, 
                       REF=SNPset[[lowBulk]]$Ref_Allele, 
                       ALT=SNPset[[lowBulk]]$ALT1, 
                       AD_REF.LOW=SNPset[[lowBulk]]$Reference,
                       AD_ALT.LOW=(SNPset[[lowBulk]]$Alternative1),
                       DP.LOW=SNPset[[lowBulk]]$RealDepth,
                       GQ.LOW=rep(NA, nrow(SNPset[[lowBulk]])) ,
                       PL.LOW=rep(NA, nrow(SNPset[[lowBulk]])),
                       SNPindex.LOW=(as.numeric(SNPset[[lowBulk]]$Alternative1))/as.numeric(SNPset[[lowBulk]]$RealDepth),
                       AD_REF.HIGH=SNPset[[highBulk]]$Reference,
                       AD_ALT.HIGH=(SNPset[[highBulk]]$Alternative1),
                       DP.HIGH=SNPset[[highBulk]]$RealDepth,
                       GQ.HIGH=rep(NA, nrow(SNPset[[highBulk]])) ,
                       PL.HIGH=rep(NA, nrow(SNPset[[highBulk]])),
                       SNPindex.HIGH=(as.numeric(SNPset[[highBulk]]$Alternative1))/as.numeric(SNPset[[highBulk]]$RealDepth), #AD_ALT.HIGH/DP.HIGH,
                       REF_FRQ=((SNPset[[highBulk]]$Reference+SNPset[[lowBulk]]$Reference)/(SNPset[[highBulk]]$RealDepth+SNPset[[lowBulk]]$RealDepth))
  )
  metaTable=dplyr::mutate(metaTable, deltaSNP=SNPindex.HIGH-SNPindex.LOW)
  return(metaTable)
}


##########################

# Now that we have a list with all comparisons for replicates (object: "replicates"), pools (object: "pools") and all 
#samples together (allPoolsInOneComparison) we can begin filtering the samples. First idea is to filter by:
#' @param SNPcomparison: lower tenth percentile > quantile(sum of depth in both bulks, na.rm = T, probs = 0.1))
#' @param maxDepthPercentile depth: higher fifth percentile > quantile(sum of depth in both bulks, na.rm = T, probs = 0.95))
#' @param minDepthPercentile sample depth - Half the min depth
#' @param windowSize
#' @param bulkSize
#' @param minimumDepthPerAllele

#This analyzer function will filter
analyzer=function(SNPcomparison, minDepthPercentile= 0.1, maxDepthPercentile=0.9, windowSize = 2.5e4, bulkSize = 20 , outlierFilt = "deltaSNP"){
  SNPcomparison=SNPcomparison %>% filter(!is.nan(deltaSNP)) %>%
    filter(!is.na(SNPindex.LOW) | !is.na(SNPindex.HIGH))%>%
    filter(!CHROM=="chrM")    #This is the last addition
  
  quants=quantile(SNPcomparison$DP.LOW + SNPcomparison$DP.HIGH, na.rm = T, probs = c(minDepthPercentile, maxDepthPercentile))
  df_filtered <- filterSNPs(
    SNPset = SNPcomparison,
    minTotalDepth = quants[1],
    maxTotalDepth = quants[2],
    #minSampleDepth = quants[1]/2,  It seems that the libraries from lung and YPD have less depth than the inoculum, so this parameter was messing up the filtering.
    verbose = TRUE)
  
  df_filtered <- runQTLseqAnalysis(df_filtered, windowSize = windowSize,
                                   popStruc = "RIL",
                                   bulkSize = bulkSize,
                                   replications = 10000, intervals = c(95, 99) )
  
  df_filtered <- runGprimeAnalysis(df_filtered, windowSize = windowSize, outlierFilter = outlierFilt  , filterThreshold = 0.1)
  
  return(df_filtered)
}


##########################


#this function will make a new table with the differences in percentages (between YPD - inoculum, and Lung - inoculum) in a sample 
#using the tabler output list as input. The sampleColumn should be the column number that shows the sample's inoculum.

percenter=function(input=allSamples, parameter="Reference", sampleColumn=1){
  #this function makes a list with the differences in percentage between the Lung/YPD samples and the inoculum. The parameter for testing is
  # either the reference or the alternative alleles.
  v=sampleColumn
  reg=input[[v]][,c(1:6)]
  if(parameter=="Reference"){reg$Lung=input[[v+1]]$Ref_percentage-input[[v]]$Ref_percentage
  reg$YPD=input[[v+2]]$Ref_percentage-input[[v]]$Ref_percentage} else
    if(parameter=="Alternative1"){reg$Lung=input[[v+1]]$Alt1_percentage-input[[v]]$Alt1_percentage
    reg$YPD=input[[v+2]]$Alt1_percentage-input[[v]]$Alt1_percentage} else
      if(parameter=="Alternative2"){reg$Lung=input[[v+1]]$Alt2_percentage-input[[v]]$Alt2_percentage
      reg$YPD=input[[v+2]]$Alt1_percentage-input[[v]]$Alt1_percentage}
  reg$CHR=as.factor(reg$CHR)
  return(reg)
}

percenter2=function(input, parameter="Reference", sampleColumn=1){
  #this function makes a list with the differences in percentage between the Lung/YPD samples and the inoculum. The parameter for testing is
  # either the reference or the alternative alleles.
  v=sampleColumn
  tester1=tabler(i=v)
  tester2=tabler(i=(v+1))
  tester3=tabler(i=(v+2))
  reg=tester1[,c(1:6)]
  if(parameter=="Reference"){reg$Lung=tester2$Ref_percentage-tester1$Ref_percentage
  reg$YPD=tester2$Ref_percentage-tester1$Ref_percentage} else
    if(parameter=="Alternative1"){reg$Lung=tester2$Alt1_percentage-tester1$Alt1_percentage
    reg$YPD=tester2$Alt1_percentage-tester1$Alt1_percentage} else
      if(parameter=="Alternative2"){reg$Lung=tester2$Alt2_percentage-tester1$Alt2_percentage
      reg$YPD=tester2$Alt1_percentage-tester1$Alt1_percentage}
  
  sampleName=colnames(tab)
  amp=str_split_fixed(sampleName, ".I", n=2)[9+v]
  
  reg$replicate=rep(amp, nrow(reg))
  rm(tester1)
  rm(tester2)
  rm(tester3)
  reg$CHR=as.factor(reg$CHR)
  reg$replicate=as.factor(reg$replicate)
  return(reg)
}

###########################
# This function will plot the enrichment of each SNP. It uses the output of the percenter() as input.

#Obs, I have to set to run the percenter() as well as define the fileNames outside the function, in the for loop to make it work.
# here what I mean:  l=percenter(sampleColumn=sampleColumn)
plotter=function(input,fileName){
  l=input
  
  ggplot(data=l) +geom_point(aes(x=POS, y=(Lung*100), colour="blue"), size=0.01)+geom_point(aes(x=POS, y=(YPD*100), colour="orange"), size=0.01)+
    scale_color_manual(name="Condition",labels = c("Lung", "YPD"), values = c("blue", "orange"))+
    facet_wrap(~factor(l$CHR, levels =unique(l$CHR) ), scales="free", ncol=3)+
    geom_hline(data=l, aes(yintercept = 0))+
    scale_y_continuous(limits = c(-100,100))+
    guides(colour = guide_legend(override.aes = list(size=6)))+
    xlab("Position")+ ylab("Enrichment (Percentage points)")+
    theme(axis.title.x=element_text(size=18), axis.title.y = element_text(size=18), axis.text.x = element_text(size=10),
          axis.text.y=element_text(size=12),strip.text.x=element_text(size = 14))
  ggsave(paste0(fileName,".png"), width = 18, units="in",height = 14)
  
  ggplot(data=l) +geom_point(aes(x=POS, y=(Lung*100), colour="blue"), size=0.001)+
    scale_color_manual(name="Condition",labels = c("Lung"), values = c("blue"))+
    facet_wrap(~factor(l$CHR, levels =unique(l$CHR) ), scales="free", ncol=3)+
    geom_hline(data=l, aes(yintercept = 0))+
    scale_y_continuous(limits = c(-100,100))+
    xlab("Position")+ ylab("Enrichment (Percentage points)")+
    theme(axis.title.x=element_text(size=18), axis.title.y = element_text(size=18), axis.text.x = element_text(size=10),
          axis.text.y=element_text(size=12),strip.text.x=element_text(size = 14))
  ggsave(paste0(fileName, ".lung_only.png"), width = 18, units="in",height = 14)
  
  ggplot(data=l) +geom_point(aes(x=POS, y=(YPD*100), colour="orange"), size=0.001)+
    scale_color_manual(name="Condition",labels = c("YPD"), values = c("Orange"))+
    facet_wrap(~factor(l$CHR, levels =unique(l$CHR) ), scales="free", ncol=3)+
    geom_hline(data=l, aes(yintercept = 0))+
    scale_y_continuous(limits = c(-100,100))+
    xlab("Position")+ ylab("Enrichment (Percentage points)")+
    theme(axis.title.x=element_text(size=18), axis.title.y = element_text(size=18), axis.text.x = element_text(size=10),
          axis.text.y=element_text(size=12),strip.text.x=element_text(size = 14))
  ggsave(paste0(fileName, ".YPD_only.png"), width = 18, units="in",height = 14)
  
}

###########################
#This is a function to make a new table with different bin sizes and the average of ∆% observed in the positions inside the bin. 
#You can play with the bin size, and it uses the output of percenter as input.  

binner=function(input, bin.size=10000){
  KN99_gen=data.frame(Chr=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chrM"),
                      size=c(2291500,1621676,1574972,1084805,1814975,1422463,1399209,1398693,1186813,1059962,1562107,774060,756017,942472,24923) )
  
  klaus=data.frame(CHR=as.character(), bin_floor=as.numeric(), bin_ceiling=as.numeric(), bin_middle=as.numeric(), 
                   lungAverage=as.numeric(), YPDaverage=as.numeric(), mut_perBin=as.numeric(), binSize=as.numeric())
  for (i in 1:nrow(KN99_gen)){
    bin.floor=0
    for (z in 1:ceiling(KN99_gen[i,2]/bin.size)){
      bin.ceiling=bin.floor+bin.size
      partial=subset(input, CHR==KN99_gen$Chr[i] & POS>bin.floor & POS<=bin.ceiling)
      lungAve=mean(partial$Lung, na.rm=T)
      YPDAve=mean(partial$YPD, na.rm=T)
      bin.middle=mean(c(bin.floor, bin.ceiling))
      klaus=rbind(klaus, c(as.character(KN99_gen$Chr[i]), as.numeric(bin.floor), as.numeric(bin.ceiling), as.numeric(bin.middle), 
                           as.numeric(lungAve), as.numeric(YPDAve), nrow(partial), bin.size),stringsAsFactors=F)
      bin.floor=bin.ceiling
    }
  } 
  colnames(klaus)=c("Chr", "binFloor","binCeiling", "binMiddle","lungMean", "YPDmean", "mut_perBin", "Bin_size")
  klaus[,2:ncol(klaus)]=sapply(klaus[,2:ncol(klaus)], as.numeric)
  return(klaus)
}

###############
binner2=function(input ,bin.size=10000){
  KN99_gen=data.frame(Chr=as.character(c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chrM")),
                      size=c(2291500,1621676,1574972,1084805,1814975,1422463,1399209,1398693,1186813,1059962,1562107,774060,756017,942472,24923) )
  input$CHROM=as.character(input$CHROM)
  
  klaus=data.frame(CHROM=as.character(), bin_floor=as.numeric(), bin_ceiling=as.numeric(), bin_middle=as.numeric(), 
                   deltaSNP=as.numeric(), triCubeDeltaSNP=as.numeric(),qvalue=as.numeric(), mut_perBin=as.numeric(), binSize=as.numeric(), 
                   Significance=as.logical())
  for (i in 1:nrow(KN99_gen)){
    bin.floor=0
    for (z in 1:ceiling(KN99_gen[i,2]/bin.size)){
      partial=data.frame()
      bin.ceiling=bin.floor+bin.size
      partial=subset(input, CHROM==KN99_gen$Chr[i] & POS>bin.floor & POS<=bin.ceiling)
      if(nrow(partial)==0){
        deltaSNP=NA
        triCubeDeltaSNP=NA
      }else if(nrow(partial)==1){
        deltaSNP=partial$deltaSNP
        triCubeDeltaSNP=partial$tricubeDeltaSNP
      }else{
        deltaSNP=mean(partial$deltaSNP, na.rm = T)
        triCubeDeltaSNP=mean(partial$tricubeDeltaSNP, na.rm = T)}
      
      significance=F
      significance=length(which(partial$qvalue<0.1))/length(partial$qvalue)>=0.5 #If half or more the snps in the window 
      if (is.na(significance)){significance=F}
      #significance=any(partial$qvalue<0.01, na.rm=T) #This is the case for window has one significant snp, all window is significant.
      
      bin.middle=mean(c(bin.floor, bin.ceiling))
      qvalue=mean(partial$qvalue)
      klaus=rbind(klaus, c(as.character(KN99_gen$Chr[i]), as.numeric(bin.floor), as.numeric(bin.ceiling), as.numeric(bin.middle), 
                           as.numeric(deltaSNP), as.numeric(triCubeDeltaSNP) ,as.numeric(qvalue) ,nrow(partial), bin.size, significance),stringsAsFactors=F)
      bin.floor=bin.ceiling
    }
  } 
  colnames(klaus)=c("CHROM", "binFloor","binCeiling", "binMiddle","deltaSNP","smoothedDeltaSNP", "qvalue", "mut_perBin", "Bin_size", "Significance")
  klaus[,2:9]=sapply(klaus[,2:9], as.numeric)
  return(klaus)
}


###############

#' This function extracts info from an imported vcf for 1 sample, outputting a nice table.
#' 
#' @param input  
#' @param bin.size The size of the window to be used in the smoothing.
#' 
#' @return a table

binner3=function(input, bin.size=10000){
  KN99_gen=data.frame(Chr=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chrM"),
                      size=c(2291500,1621676,1574972,1084805,1814975,1422463,1399209,1398693,1186813,1059962,1562107,774060,756017,942472,24923) )
  
  smoothedDeltaSNP=vector()
  for (i in 1:nrow(KN99_gen)){
    bin.floor=0
    for (z in 1:ceiling(KN99_gen[i,2]/bin.size)){
      bin.ceiling=bin.floor+bin.size
      partial=subset(input, CHROM==KN99_gen$Chr[i] & POS>bin.floor & POS<=bin.ceiling)
      tempValue=vector()
     if(nrow(partial)==1){tempValue=partial$deltaSNP}
      else if (nrow(partial)>1){tempValue=tricubeStat(POS=partial$POS, Stat=partial$deltaSNP, windowSize = bin.size)}
      if(length(tempValue)>0){smoothedDeltaSNP=c(smoothedDeltaSNP, tempValue) }   
      bin.floor=bin.ceiling
    }
  } 

  return(smoothedDeltaSNP)
}


###############

#' Calculate tricube weighted statistics for each SNP
#'
#' Uses local regression (wrapper for \code{\link[locfit]{locfit}}) to predict a
#' tricube smoothed version of the statistic supplied for each SNP. This works as a
#' weighted average across neighboring SNPs that accounts for Linkage
#' disequilibrium (LD) while minizing noise attributed to SNP calling errors.
#' Values for neighboring SNPs within the window are weighted by physical
#' distance from the focal SNP.
#'
#' @return Returns a vector of the weighted statistic caluculted with a
#'   tricube smoothing kernel
#'
#' @param POS A vector of genomic positions for each SNP
#' @param Stat A vector of values for a given statistic for each SNP
#' @param WinSize the window size (in base pairs) bracketing each SNP for which
#'   to calculate the statitics. Magwene et. al recommend a window size of ~25
#'   cM, but also recommend optionally trying several window sizes to test if
#'   peaks are over- or undersmoothed.
#' @examples df_filt_4mb$Gprime <- tricubeStat(POS, Stat = GStat, WinSize = 4e6)
#' @seealso \code{\link{getG}} for G statistic calculation
#' @seealso \code{\link[locfit]{locfit}} for local regression

tricubeStat <- function(POS, Stat, windowSize = 2e6)
{
  if (windowSize <= 0)
    stop("A positive smoothing window is required")
  stats::predict(locfit::locfit(Stat ~ locfit::lp(POS, h = windowSize, deg = 0)), POS)
}

#################


#This function makes the distribution plots for the bins, using the output table of the binner()palette() function as input.  

binDistPlotter=function(inputBinTable, windowsSizes=c(5000,10000, 20000, 50000), fileName ){
  binnerDataframe=data.frame(CHR=as.character(), bin_floor=as.numeric(), bin_ceiling=as.numeric(), bin_middle=as.numeric(), 
                             lungAverage=as.numeric(), YPDaverage=as.numeric(), mut_perBin=as.numeric(), bin.size=as.numeric())
  for (z in 1:length(windowsSizes)){
    binnerDataframe=rbind(binnerDataframe, binner(inputBinTable, bin.size = windowsSizes[z]))
  }
  zeta=melt(binnerDataframe, id.vars = colnames(binnerDataframe[,c(1:4,7,8)]))
  zeta_lung=subset(zeta,variable=="lungMean")
  zeta_ypd=subset(zeta,variable=="YPDmean")
  
  ggplot(data=zeta_lung, aes(x=(zeta_lung$value*100), color=factor(zeta_lung$Bin_size)))+geom_density()+
    #facet_wrap(~factor(zeta_lung$Chr,levels=unique(zeta_lung$Chr)), ncol = 3)+
    scale_x_continuous(limits = c(-100,100))+
    scale_color_manual(name="Bin size",labels = as.character(windowsSizes) , values = palette()[1:length(windowsSizes)])+
    theme(axis.title.x=element_text(size=12), axis.title.y = element_text(size=12), axis.text.x = element_text(size=8),
          axis.text.y=element_text(size=6),strip.text.x=element_text(size = 10))+
    xlab("Enrichment (Percentage points)")
  ggsave(paste0(fileName,".binLungDistribution.png"))
  
  ggplot(data=zeta_ypd, aes(x=(zeta_ypd$value*100), color=factor(zeta_ypd$Bin_size)))+geom_density()+
    #facet_wrap(~factor(zeta_lung$Chr,levels=unique(zeta_lung$Chr)), ncol = 3)+
    scale_x_continuous(limits = c(-100,100))+
    scale_color_manual(name="Bin size",labels = as.character(windowsSizes) , values = palette()[1:length(windowsSizes)])+
    theme(axis.title.x=element_text(size=12), axis.title.y = element_text(size=12), axis.text.x = element_text(size=8),
          axis.text.y=element_text(size=6),strip.text.x=element_text(size = 10))+
    xlab("Enrichment (Percentage points)")
  ggsave(paste0(fileName,".binYpdDistribution.png"))
}

###############

# Function to set comparisons between two samples and output the table in a way that the QTLseq package can use.

#' Uses two element of a list to create a table that can be used by the package QTLseq.
#' 
#' @param SNPset A list object. Each eelement of the list is a table showing data from one sample, the table is set as the output of the tabler function.
#' @param lowBulk String. The name of the sample (and of the element in the SNPset list) corresponding to the Low Bulk.
#' @param highBulk String. The name of the sample (and of the element in the SNPset list) corresponding to the High Bulk.
#' 
#' @return a data.frame with the format that can be understood by QTLseq. 

picker=function(SNPset, lowBulk, highBulk){
  if(!lowBulk %in% names(SNPset) | !highBulk %in% names(SNPset)){
    stop("Confirm if the bulks' names are present in the SNP set (vcf file) provided.")
  }
  metaTable=data.frame(CHROM=SNPset[[lowBulk]]$CHR, POS=SNPset[[lowBulk]]$POS, REF=SNPset[[lowBulk]]$Ref_allele, 
                       ALT=SNPset[[lowBulk]]$ALT1, 
                       AD_REF.LOW=SNPset[[lowBulk]]$Reference,
                       AD_ALT.LOW=(SNPset[[lowBulk]]$Alternative1+(SNPset[[lowBulk]]$Alternative2)),
                       DP.LOW=SNPset[[lowBulk]]$realDepth,
                       GQ.LOW=rep(NA, nrow(SNPset[[lowBulk]])) ,
                       PL.LOW=rep(NA, nrow(SNPset[[lowBulk]])),
                       SNPindex.LOW=(SNPset[[lowBulk]]$Alternative1+(SNPset[[lowBulk]]$Alternative2))/SNPset[[lowBulk]]$realDepth,
                       AD_REF.HIGH=SNPset[[highBulk]]$Reference,
                       AD_ALT.HIGH=(SNPset[[highBulk]]$Alternative1+(SNPset[[highBulk]]$Alternative2)),
                       DP.HIGH=SNPset[[highBulk]]$realDepth,
                       GQ.HIGH=rep(NA, nrow(SNPset[[highBulk]])) ,
                       PL.HIGH=rep(NA, nrow(SNPset[[highBulk]])),
                       SNPindex.HIGH=(SNPset[[highBulk]]$Alternative1+(SNPset[[highBulk]]$Alternative2))/SNPset[[highBulk]]$realDepth, #AD_ALT.HIGH/DP.HIGH,
                       REF_FRQ=((SNPset[[highBulk]]$Reference+SNPset[[lowBulk]]$Reference)/(SNPset[[highBulk]]$realDepth+SNPset[[lowBulk]]$realDepth)))
  
  metaTable=dplyr::mutate(metaTable, deltaSNP=SNPindex.HIGH-SNPindex.LOW)
  
}


##############################


# Now that we have a list with all comparisons for replicates (object: "replicates"), pools (object: "pools") and all 
#samples together (allPoolsInOneComparison) we can begin filtering the samples. First idea is to filter by:
#' @param MinDepthPercentile: lower tenth percentile > quantile(sum of depth in both bulks, na.rm = T, probs = 0.1))
#' @param Max depth: highr fith percentile > quantile(sum of depth in both bulks, na.rm = T, probs = 0.95))
#' @param Min sample depth - Half the min depth

#This analyzer function will filter
analyzer=function(SNPcomparison, minDepthPercentile= 0.1, maxDepthPercentile=0.9, windowSize = 2.5e4, bulkSize = 20 , outlierFilt = "deltaSNP"){
  
  quants=quantile(SNPcomparison$DP.LOW + SNPcomparison$DP.HIGH, na.rm = T, probs = c(minDepthPercentile, maxDepthPercentile))
  df_filtered <- filterSNPs(
    SNPset = SNPcomparison,
    minTotalDepth = quants[1],
    maxTotalDepth = quants[2],
    #minSampleDepth = quants[1]/2,  It seems that the libraries from lung and YPD have less depth than the inoculum, so this parameter was messing up the filtering.
    verbose = TRUE)
  
  df_filtered <- runQTLseqAnalysis(df_filtered, windowSize = windowSize,
                                   popStruc = "RIL",
                                   bulkSize = bulkSize,
                                   replications = 10000, intervals = c(95, 99) )
  
  df_filtered <- runGprimeAnalysis(df_filtered, windowSize = windowSize, outlierFilter = outlierFilt  , filterThreshold = 0.1)
  
  return(df_filtered)
}


######################################################

simulateConfInt <-
  function(SNPset, popStruc = "F2",
           bulkSize,
           depth = 1:100,
           replications = 10000,
           filter = 0.3,
           intervals = c(0.05, 0.025)) {
    if (popStruc == "F2") {
      message(
        "Assuming bulks selected from F2 population, with ",
        as.character(bulkSize),
        " individuals per bulk."
      )
    } else {
      message(
        "Assuming bulks selected from RIL population, with ",
        as.character(bulkSize),
        " individuals per bulk."
      )
    }
    
    #makes a vector of possible alt allele frequencies once. this is then sampled for each replicate
    
    #This part commented below was replaced in the call below.
    
    #tmp_freq <-
     # replicate(n = replications * 10, simulateAlleleFreq(n = bulkSize, pop = popStruc))
    
    message(
      paste0(
        "Simulating ",
        as.character(replications),
        " SNPs with reads at each depth: ",
        as.character(min(depth)),
        "-",
        as.character(max(depth))
      )
    )
    message(paste0(
      "Keeping SNPs with >= ",
      as.character(filter),
      " SNP-index in both simulated bulks"
    ))
    
    # tmp allele freqs are sampled to produce 'replicate' numbers of probablities. these 
    # are then used as altFreq probs to simulate SNP index values, per bulk.
    CI <- sapply(
      X = depth,
      FUN = function(x)
      {
        quantile(
          x = simulateSNPindex(
            depth = x,
            altFreq1 = SNPset$SNPindex.LOW, #I adapted this. Instead of sampling the tmp_freq, we are using the SNP index from the inoculum to calculate this.
            
            altFreq2 =  SNPset$SNPindex.LOW,
            replicates = replications,
            filter = filter
          ),
          probs = intervals,
          names = TRUE
        )
      }
    )
    
    CI <- as.data.frame(CI)
    
    if (length(CI) > 1) {
      CI <- data.frame(t(CI))
    }
    
    names(CI) <- paste0("CI_", 100 - (intervals * 200))
    CI <- cbind(depth, CI)
    
    #to long format for easy plotting
    # tidyr::gather(data = CI,
    #     key = interval,
    #     convert = TRUE,
    #     value = SNPindex,-depth) %>%
    #     dplyr::mutate(Confidence = factor(ifelse(
    #         interval > 0.5,
    #         paste0(round((1 - interval) * 200, digits = 1), "%"),
    #         paste0((interval * 200), "%")
    # )))
    CI
  }


######################################################


runQTLseqAnalysis <- function(SNPset, windowSize = 1e6,
                              popStruc = "F2",
                              bulkSize,
                              depth = NULL,
                              replications = 10000,
                              filter = 0.3,
                              intervals = c(95, 99)
) {
  
  message("Counting SNPs in each window...")
  SNPset <- SNPset %>%
    dplyr::group_by(CHROM) %>%
    dplyr::mutate(nSNPs = countSNPs_cpp(POS = POS, windowSize = windowSize))
  
  message("Calculating tricube smoothed delta SNP index...")
  SNPset <- SNPset %>%
    dplyr::mutate(tricubeDeltaSNP = tricubeStat(POS = POS, Stat = deltaSNP, windowSize))
  
  #convert intervals to quantiles
  if (all(intervals >= 1)) {
    message("Returning the following two sided confidence intervals: ", paste(intervals, collapse = ", "))
    quantiles <- (100 - intervals) / 200
  } else {
    stop(
      "Convidence intervals ('intervals' paramater) should be supplied as two-sided percentiles. i.e. If intervals = '95' will return the two sided 95% confidence interval, 2.5% on each side."
    )
  }
  
  #calculate min depth per snp between bulks
  SNPset <-
    SNPset %>% 
    dplyr::mutate(minDP = pmin(DP.LOW, DP.HIGH))
  
  SNPset <-
    SNPset %>% 
    dplyr::group_by(CHROM) %>% 
    dplyr::mutate(tricubeDP = floor(tricubeStat(POS, minDP, windowSize = windowSize)))
  
  if (is.null(depth)) {
    message(
      "Variable 'depth' not defined, using min and max depth from data: ",
      min(SNPset$minDP),
      "-",
      max(SNPset$minDP)
    )
    depth <- min(SNPset$minDP):max(SNPset$minDP)
  }
  
  #simulate confidence intervals
  CI <-
    simulateConfInt(SNPset,
      popStruc = popStruc,
      bulkSize = bulkSize,
      depth = depth,
      replications = replications,
      filter = filter,
      intervals = quantiles
    )
  
  
  #match name of column for easier joining of repeat columns
  names(CI)[1] <- "tricubeDP"
  
  #use join as a quick way to match min depth to matching conf intervals.
  SNPset <-
    dplyr::left_join(x = SNPset,
                     y = CI #, commented out because of above change. need to remove eventually
                     # by = c("tricubeDP" = "depth")
    )
  
  as.data.frame(SNPset)
  
}


######################

# Now that we have a list with all comparisons for replicates (object: "replicates"), pools (object: "pools") and all 
#samples together (allPoolsInOneComparison) we can begin filtering the samples. First idea is to filter by:
#' @param MinDepthPercentile: lower tenth percentile > quantile(sum of depth in both bulks, na.rm = T, probs = 0.1))
#' @param Max depth: highr fith percentile > quantile(sum of depth in both bulks, na.rm = T, probs = 0.95))
#' @param Min sample depth - Half the min depth

#This analyzer function will filter
analyzer=function(SNPcomparison, minDepthPercentile= 0.1, maxDepthPercentile=0.9, windowSize = 2.5e4, bulkSize = 20 , outlierFilt = "deltaSNP"){
  
  SNPcomparison=SNPcomparison %>% filter(!is.nan(deltaSNP)) %>%
    filter(!is.na(SNPindex.LOW) | !is.na(SNPindex.HIGH))  #This is the last addition
  
  quants=quantile(SNPcomparison$DP.LOW + SNPcomparison$DP.HIGH, na.rm = T, probs = c(minDepthPercentile, maxDepthPercentile))
  df_filtered <- filterSNPs(
    SNPset = SNPcomparison,
    minTotalDepth = quants[1],
    maxTotalDepth = quants[2],
    #minSampleDepth = quants[1]/2,  It seems that the libraries from lung and YPD have less depth than the inoculum, so this parameter was messing up the filtering.
    verbose = TRUE)
  
  df_filtered <- runQTLseqAnalysis(df_filtered, windowSize = windowSize,
                                   popStruc = "RIL",
                                   bulkSize = bulkSize,
                                   replications = 10000, intervals = c(95, 99) )
  
  df_filtered <- runGprimeAnalysis(df_filtered, windowSize = windowSize, outlierFilter = outlierFilt  , filterThreshold = 0.1)
  
  return(df_filtered)
}

###################################

#' Simulates a delta SNP-index with replication
#'
#' @param depth integer. A read depth for which to replicate SNP-index calls.
#' @param altFreq1 numeric. The alternate allele frequency for bulk A. 
#' @param altFreq2 numeric. The alternate allele frequency for bulk B. 
#' @param replicates integer. The number of bootstrap replications.
#' @param filter numeric. an optional minimum SNP-index filter
#'
#' @return Returns a vector of length replicates delta SNP-indeces 


simulateSNPindex <-
  function(depth,
           altFreq1,
           altFreq2,
           replicates = 10000,
           filter = NULL) {
    
    SNPindex_H <- rbinom(replicates, size = depth, altFreq1) / depth
    SNPindex_L <- rbinom(replicates, size = depth, altFreq2) / depth
    deltaSNP <- SNPindex_H - SNPindex_L
    
    if (!is.null(filter)) {
      deltaSNP <- deltaSNP[SNPindex_H >= filter | SNPindex_L >= filter]
    }
    deltaSNP
  }

##########################################################################################################################################

#Functions for calculating and manipulating the G statistic

#' Calculates the G statistic
#'
#' The function is used by \code{\link{runGprimeAnalysis}} to calculate the G
#' statisic G is defined by the equation: \deqn{G = 2*\sum_{i=1}^{4}
#' n_{i}*ln\frac{obs(n_i)}{exp(n_i)}}{G = 2 * \sum n_i * ln(obs(n_i)/exp(n_i))}
#' Where for each SNP, \eqn{n_i} from i = 1 to 4 corresponds to the reference
#' and alternate allele depths for each bulk, as described in the following
#' table: \tabular{rcc}{ Allele \tab High Bulk \tab Low Bulk \cr Reference \tab
#' \eqn{n_1} \tab \eqn{n_2} \cr Alternate \tab \eqn{n_3} \tab \eqn{n_4} \cr}
#' ...and \eqn{obs(n_i)} are the observed allele depths as described in the data
#' frame. Method 1 calculates the G statistic using expected values assuming
#' read depth is equal for all alleles in both bulks: \deqn{exp(n_1) = ((n_1 +
#' n_2)*(n_1 + n_3))/(n_1 + n_2 + n_3 + n_4)} \deqn{exp(n_2) = ((n_2 + n_1)*(n_2
#' + n_4))/(n_1 + n_2 + n_3 + n_4)} etc...
#'
#' @param LowRef A vector of the reference allele depth in the low bulk
#' @param HighRef A vector of the reference allele depth in the high bulk
#' @param LowAlt A vector of the alternate allele depth in the low bulk
#' @param HighAlt A vector of the alternate allele depth in the high bulk
#'
#' @return A vector of G statistic values with the same length as
#'
#' @seealso \href{https://doi.org/10.1371/journal.pcbi.1002255}{The Statistics
#'   of Bulk Segregant Analysis Using Next Generation Sequencing}
#'   \code{\link{tricubeStat}} for G prime calculation

getG <- function(LowRef, HighRef, LowAlt, HighAlt)
{
  exp <- c(
    (LowRef + HighRef) * (LowRef + LowAlt) / (LowRef + HighRef + LowAlt + HighAlt),
    (LowRef + HighRef) * (HighRef + HighAlt) / (LowRef + HighRef + LowAlt + HighAlt),
    (LowRef + LowAlt) * (LowAlt + HighAlt) / (LowRef + HighRef + LowAlt + HighAlt),
    (LowAlt + HighAlt) * (HighRef + HighAlt) / (LowRef + HighRef + LowAlt + HighAlt)
  )
  obs <- c(LowRef, HighRef, LowAlt, HighAlt)
  
  G <-
    2 * (rowSums(obs * log(
      matrix(obs, ncol = 4) / matrix(exp, ncol = 4)
    )))
  return(G)
}
##########################################################################################################################################


#' Identifies which variants are statistically significant in x number of pools.
#'
#' @param depth integer. A read depth for which to replicate SNP-index calls.
#' @param altFreq1 numeric. The alternate allele frequency for bulk A. 
#' @param altFreq2 numeric. The alternate allele frequency for bulk B. 
#' @param replicates integer. The number of bootstrap replications.
#' @param filter numeric. an optional minimum SNP-index filter
#'
#' @return Returns a vector of length replicates delta SNP-indeces 
#' 
plotMakerBSA=function(datum, xplot=datum$binMiddle/1000, yplot=datum$All_pools){
  ggplot(data=datum,aes(x=xplot, y=yplot))+
    geom_pointrange(data=temp,
                    aes(ymin=pool_min, 
                        ymax=pool_max ,
                        x=binMiddle/1000, 
                        y=All_pools,
                        size=combine,
                        colour=combine,
                        alpha=combine, 
                        order="combine"))+
    scale_size_manual(name="Condition and significance",
                      labels=c("Lung, not significant","YPD, not significant","Lung, significant",  "YPD, significant"),
                      values = c(0.1,0.1,0.4, 0.4), drop=F)+
    scale_colour_manual(name="Condition and significance",
                        labels=c("Lung, not significant","YPD, not significant","Lung, significant",  "YPD, significant"),
                        values = c('FALSE_Lungs'="red",'FALSE_YPD'="black",
                                   'TRUE_Lungs'="red",'TRUE_YPD'="black"), drop=F)+
    scale_alpha_manual(name="Condition and significance",
                       labels=c("Lung, not significant","YPD, not significant","Lung, significant",  "YPD, significant"),
                       values = c(0.95,0.95,1, 1), drop=F)+
    #geom_vline(xintercept  = IR_list[[paste0("chr",i)]],color="black", size=0.3, linetype="dashed")+
    #geom_vline(xintercept = 0, size=0.3)+
    labs(title=paste0("Chromosome ", i), x=paste0("Position on chromosome ",i," (kb)"),y="Change in allele frequency (BSA)", color="Pool")+
    geom_hline(yintercept = c(0), color="black", size=0.6)+
    scale_x_continuous(breaks = seq(from=0, to=2400, by=300), limits=c(0, 2400),expand = c(0, 0))+
    scale_y_continuous(breaks = c(-0.5, -0.25, 0,0.25,0.5), limits = c(-0.55,0.55), expand = c(0,0))+
    theme(#legend.title = element_text(size=18), 
      legend.position = "right", 
      plot.title = element_blank(),
      legend.text = element_text(size=16),
      axis.title=element_text(size=32),
      axis.title.y =element_text(margin = margin(r=25, t=0, l=5, b=0)),
      axis.title.x =element_text(margin = margin(r=0, t=25, l=0, b=5)),
      axis.text=element_text(size = 28),
      axis.ticks.length = unit(0.3,"cm"),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      panel.background = element_rect(fill = "white",colour = "black", size=1),
      text=element_text(family="sans"))
  
  
  
}


#' Non-parametric estimation of the null distribution of G'
#'
#' The function is used by \code{\link{runGprimeAnalysis}} to estimate p-values for the
#' weighted G' statistic based on the non-parametric estimation method described
#' in Magwene et al. 2011. Breifly, using the natural log of Gprime a median
#' absolute deviation (MAD) is calculated. The Gprime set is trimmed to exclude
#' outlier regions (i.e. QTL) based on Hampel's rule. An alternate method for
#' filtering out QTL is proposed using absolute delta SNP indeces greater than
#' a set threshold to filter out potential QTL. An estimation of the mode of the trimmed set
#' is calculated using the \code{\link[modeest]{mlv}} function from the package
#' modeest. Finally, the mean and variance of the set are estimated using the
#' median and mode and p-values are estimated from a log normal distribution.
#'
#' @param Gprime a vector of G prime values (tricube weighted G statistics)
#' @param deltaSNP a vector of delta SNP values for use for QTL region filtering
#' @param outlierFilter one of either "deltaSNP" or "Hampel". Method for
#'   filtering outlier (ie QTL) regions for p-value estimation
#' @param filterThreshold The absolute delta SNP index to use to filter out putative QTL
#' @export getPvals

getPvals <-
  function(Gprime,
           deltaSNP = NULL,
           outlierFilter = c("deltaSNP", "Hampel"),
           filterThreshold)
  {
    
    if (outlierFilter == "deltaSNP") {
      
      if (abs(filterThreshold) >= 0.5) {
        stop("filterThreshold should be less than 0.5")
      }
      
      message("Using deltaSNP-index to filter outlier regions with a threshold of ", filterThreshold)
      trimGprime <- Gprime[abs(deltaSNP) < abs(filterThreshold)]
    } else {
      message("Using Hampel's rule to filter outlier regions")
      lnGprime <- log(Gprime)
      
      medianLogGprime <- median(lnGprime)
      
      # calculate left median absolute deviation for the trimmed G' prime set
      MAD <-
        median(medianLogGprime - lnGprime[lnGprime <= medianLogGprime])
      
      # Trim the G prime set to exclude outlier regions (i.e. QTL) using Hampel's rule
      trimGprime <-
        Gprime[lnGprime - median(lnGprime) <= 5.2 * MAD]
    }
    
    medianTrimGprime <- median(trimGprime)
    
    # estimate the mode of the trimmed G' prime set using the half-sample method
    message("Estimating the mode of a trimmed G prime set using the 'modeest' package...")
    modeTrimGprime <-
      modeest::mlv(x = trimGprime, bw = 0.5, method = "hsm")[1]
    
    muE <- log(medianTrimGprime)
    varE <- abs(muE - log(modeTrimGprime))
    #use the log normal distribution to get pvals
    message("Calculating p-values...")
    pval <-
      1 - plnorm(q = Gprime,
                 meanlog = muE,
                 sdlog = sqrt(varE))
    
    return(pval)
  }


#########################################
#' Identify QTL using a smoothed G statistic
#'
#' A wrapper for all the functions that perform the full G prime analysis to
#' identify QTL. The following steps are performed:\cr 1) Genome-wide G
#' statistics are calculated by \code{\link{getG}}. \cr G is defined by the
#' equation: \deqn{G = 2*\sum_{i=1}^{4} n_{i}*ln\frac{obs(n_i)}{exp(n_i)}}{G = 2
#' * \sum n_i * ln(obs(n_i)/exp(n_i))} Where for each SNP, \eqn{n_i} from i = 1
#' to 4 corresponds to the reference and alternate allele depths for each bulk,
#' as described in the following table: \tabular{rcc}{ Allele \tab High Bulk
#' \tab Low Bulk \cr Reference \tab \eqn{n_1} \tab \eqn{n_2} \cr Alternate \tab
#' \eqn{n_3} \tab \eqn{n_4} \cr} ...and \eqn{obs(n_i)} are the observed allele
#' depths as described in the data frame. \code{\link{getG}} calculates the G statistic
#' using expected values assuming read depth is equal for all alleles in both
#' bulks: \deqn{exp(n_1) = ((n_1 + n_2)*(n_1 + n_3))/(n_1 + n_2 + n_3 + n_4)}
#' \deqn{exp(n_2) = ((n_2 + n_1)*(n_2 + n_4))/(n_1 + n_2 + n_3 + n_4)}
#' \deqn{exp(n_3) = ((n_3 + n_1)*(n_3 + n_4))/(n_1 + n_2 + n_3 + n_4)}
#' \deqn{exp(n_4) = ((n_4 + n_2)*(n_4 + n_3))/(n_1 + n_2 + n_3 + n_4)}\cr 2) G'
#' - A tricube-smoothed G statistic is predicted by local regression within each
#' chromosome using \code{\link{tricubeStat}}. This works as a weighted average
#' across neighboring SNPs that accounts for Linkage disequilibrium (LD) while
#' minizing noise attributed to SNP calling errors. G values for neighboring
#' SNPs within the window are weighted by physical distance from the focal SNP.
#' \cr \cr 3) P-values are estimated based using the non-parametric method
#' described by Magwene et al. 2011 with the function \code{\link{getPvals}}.
#' Breifly, using the natural log of Gprime a median absolute deviation (MAD) is
#' calculated. The Gprime set is trimmed to exclude outlier regions (i.e. QTL)
#' based on Hampel's rule. An alternate method for filtering out QTL is proposed
#' using absolute delta SNP indeces greater than 0.1 to filter out potential
#' QTL. An estimation of the mode of the trimmed set is calculated using the
#' \code{\link[modeest]{mlv}} function from the package modeest. Finally, the
#' mean and variance of the set are estimated using the median and mode and
#' p-values are estimated from a log normal distribution. \cr \cr 4) Negative
#' Log10- and Benjamini-Hochberg adjusted p-values are calculated using
#' \code{\link[stats]{p.adjust}}
#'
#' @param SNPset Data frame SNP set containing previously filtered SNPs
#' @param windowSize the window size (in base pairs) bracketing each SNP for which
#'   to calculate the statitics. Magwene et. al recommend a window size of ~25
#'   cM, but also recommend optionally trying several window sizes to test if
#'   peaks are over- or undersmoothed.
#' @param outlierFilter one of either "deltaSNP" or "Hampel". Method for
#'   filtering outlier (ie QTL) regions for p-value estimation
#' @param filterThreshold The absolute delta SNP index to use to filter out putative QTL (default = 0.1)
#' @param ... Other arguments passed to \code{\link[locfit]{locfit}} and
#'   subsequently to \code{\link[locfit]{locfit.raw}}() (or the lfproc). Usefull
#'   in cases where you get "out of vertex space warnings"; Set the maxk higher
#'   than the default 100. See \code{\link[locfit]{locfit.raw}}(). But if you
#'   are getting that warning you should seriously consider increasing your
#'   window size.
#'   
#' @return The supplied SNP set tibble after G' analysis. Includes five new
#'   columns: \itemize{\item{G - The G statistic for each SNP} \item{Gprime -
#'   The tricube smoothed G statistic based on the supplied window size}
#'   \item{pvalue - the pvalue at each SNP calculatd by non-parametric
#'   estimation} \item{negLog10Pval - the -Log10(pvalue) supplied for quick
#'   plotting} \item{qvalue - the Benajamini-Hochberg adjusted p-value}}
#'
#'
#' @importFrom dplyr %>%
#'
#' @export runGprimeAnalysis
#'
#' @examples df_filt <- runGprimeAnalysis(df_filt,windowSize = 2e6,outlierFilter = "deltaSNP")
#' @useDynLib QTLseqr
#' @importFrom Rcpp sourceCpp


runGprimeAnalysis <-
  function(SNPset,
           windowSize = 1e6,
           outlierFilter = "deltaSNP",
           filterThreshold = 0.1, 
           ...)
  {
    message("Counting SNPs in each window...")
    SNPset <- SNPset %>%
      dplyr::group_by(CHROM) %>%
      dplyr::mutate(nSNPs = countSNPs_cpp(POS = POS, windowSize = windowSize))
    
    message("Calculating tricube smoothed delta SNP index...")
    SNPset <- SNPset %>%
      dplyr::mutate(tricubeDeltaSNP = tricubeStat(POS = POS, Stat = deltaSNP, windowSize, ...))
    
    message("Calculating G and G' statistics...")
    SNPset <- SNPset %>%
      dplyr::mutate(
        G = getG(
          LowRef = AD_REF.LOW,
          HighRef = AD_REF.HIGH,
          LowAlt = AD_ALT.LOW,
          HighAlt = AD_ALT.HIGH
        ),
        Gprime = tricubeStat(
          POS = POS,
          Stat = G,
          windowSize = windowSize,
          ...
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        pvalue = getPvals(
          Gprime = Gprime,
          deltaSNP = deltaSNP,
          outlierFilter = outlierFilter,
          filterThreshold = filterThreshold
        ),
        negLog10Pval = -log10(pvalue),
        qvalue = p.adjust(p = pvalue, method = "BH")
      )
    
    return(as.data.frame(SNPset))
  }
