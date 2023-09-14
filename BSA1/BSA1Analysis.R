library(tidyverse)

source("BSA6tools.R")

#####Parent strains -
parents=list()
for(x in c("KN99a", "TDY1993")){
  parents[[x]]=read.table(file=paste0("VCFtables/",x,".txt"),header = T,stringsAsFactors = F, numerals = "allow.loss");
  parents[[x]]$Alternative1=as.integer(parents[[x]]$Alternative1)
}
#Finding the positions where KN99a is reference and TDY1993 (C8) is alternative.
parents=parents[[1]][parents[[2]]$Filt_Genotype=="Alternative" & parents[[1]]$Filt_Genotype=="Reference",] 

poolsBSA2=list.files(path="Parents_tables") %>%
  sub(pattern=".txt", replacement = "")

mySamples=list()
for (i in poolsBSA2){
  mySamples[[i]]=read.table(file = paste0("VCF_tables/",i,
                                          ".txt"), header = T,stringsAsFactors = F, numerals = "allow.loss");
  mySamples[[i]]=mySamples[[i]][paste(mySamples[[i]]$CHR, mySamples[[i]]$POS, sep="_") %in% paste(parents$CHR, parents$POS, sep="_") ,]
  mySamples[[i]]$Alternative1=as.integer(mySamples[[i]]$Alternative1)
}

#Summary of data from the coverage
SampleMean=sapply(mySamples, FUN=function(x){mean(x$RealDepth)})
SampleMedian=sapply(mySamples, FUN=function(x){median(x$RealDepth)})

pools=c(paste0("P",c(1:5),".1"),paste0("P",c(1:5),".2")) #The different pools, in this case, P1 to P5 and all
Inoculum=grep(pattern="P[[:alnum:]].{1,3}I", poolsBSA2, value=T) #The Inoculum samples
YPD=grep(pattern="P[[:alnum:]].{1,3}Y", poolsBSA2, value=T) # The YPD samples
Lungs=grep(pattern="P[[:alnum:]].{1,3}L", poolsBSA2, value=T) # The Lung samples
Brains=grep(pattern="P[[:alnum:]].{1,3}B", poolsBSA2, value=T) # The Brain samples

replicatesBSA=list()
for(i in pools){
  for(z in grep(pattern = i, x=c(YPD,Lungs,Brains), value = T)){
    lowbulk=grep(pattern = i, x = Inoculum, value = T)
    highbulk=z
    replicatesBSA[[highbulk]]=picker2(mySamples, lowBulk = lowbulk, highBulk = highbulk)
    print(paste("Comparisons made:",highbulk,"x", lowbulk, sep = " "))}
}


#Making a new list with the sums of all replicates of a pool

summedPoolsBSA=list() # a list with the summed depth and allele counts for all replicates of a specific pool
poolsBSA=list() # a list of comparisons of pooled replicates (all the replicates from a specific pool) for each condition using picker()
pools=paste0("P",1:5)

for (i in pools) {
  #Inoculum
  meta=mySamples[[grep(pattern = i, x=Inoculum, value = T)[1]]][,1:4] #this is the basis for the table; the four first columns of mySamples.
  meta$RealDepth=rowSums(sapply(grep(pattern = i, x=Inoculum, value = T) ,FUN = function(x){cbind(mySamples[[x]][,8])}))
  meta$Reference=rowSums(sapply(grep(pattern = i, x=Inoculum, value = T) ,FUN = function(x){cbind(mySamples[[x]][,9])}))
  meta$Alternative1=rowSums(sapply(grep(pattern = i, x=Inoculum, value =T) ,FUN = function(x){cbind(mySamples[[x]][,10])}))
  summedPoolsBSA[[paste0(i,"I")]]=meta
  rm(meta)
  
  #YPD
  meta=mySamples[[grep(pattern = i, x=YPD, value = T)[1]]][,1:4] #this is the basis for the table; the four first columns of mySamples.
  meta$RealDepth=rowSums(sapply(grep(pattern = i, x=YPD, value = T) ,FUN = function(x){cbind(mySamples[[x]][,8])}))
  meta$Reference=rowSums(sapply(grep(pattern = i, x=YPD, value = T) ,FUN = function(x){cbind(mySamples[[x]][,9])}))
  meta$Alternative1=rowSums(sapply(grep(pattern = i, x=YPD, value = T) ,FUN = function(x){cbind(mySamples[[x]][,10])}))
  summedPoolsBSA[[paste0(i,"Y")]]=meta
  rm(meta)
  
  #Lungs
  meta=mySamples[[grep(pattern = i, x=Lungs, value = T)[1]]][,1:4] #this is the basis for the table; the four first columns of mySamples.
  meta$RealDepth=rowSums(sapply(grep(pattern = i, x=Lungs, value = T) ,FUN = function(x){cbind(mySamples[[x]][,8])}))
  meta$Reference=rowSums(sapply(grep(pattern = i, x=Lungs, value = T) ,FUN = function(x){cbind(mySamples[[x]][,9])}))
  meta$Alternative1=rowSums(sapply(grep(pattern = i, x=Lungs, value = T) ,FUN = function(x){cbind(mySamples[[x]][,10])}))
  summedPoolsBSA[[paste0(i,"L")]]=meta
  rm(meta)
  
  meta=mySamples[[grep(pattern = i, x=Brains, value = T)[1]]][,1:4] #this is the basis for the table; the four first columns of mySamples.
  meta$RealDepth=rowSums(sapply(grep(pattern = i, x=Brains, value = T) ,FUN = function(x){cbind(mySamples[[x]][,8])}))
  meta$Reference=rowSums(sapply(grep(pattern = i, x=Brains, value = T) ,FUN = function(x){cbind(mySamples[[x]][,9])}))
  meta$Alternative1=rowSums(sapply(grep(pattern = i, x=Brains, value = T) ,FUN = function(x){cbind(mySamples[[x]][,10])}))
  summedPoolsBSA[[paste0(i,"B")]]=meta
  rm(meta)
}

##Making the contrasts of low and high bulk with the summed pools

for(i in pools){ #The pools
  for(y in c("Y","L","B")){
    lowbulk=paste0(i,"I")
    highbulk=paste0(i,y)
    if(lowbulk %in% names(summedPoolsBSA) & highbulk %in% names(summedPoolsBSA) ){
      poolsBSA[[highbulk]]=picker2(summedPoolsBSA, lowBulk = lowbulk, highBulk = highbulk)
      print(paste("Comparisons made:",highbulk,"x", lowbulk, sep = " "))}
  }
}

#Here we start summing all pools together for each condition, and then making the comparisons.

allPoolsInOneBSA=list() #This will be the list containing three itens, the sum of all inocula, all YPDs and all Inoculum.
allPoolsInOneComparisonBSA=list()

for (y in c("I","Y","L", "B")){ #This loop goes around for each condition.
  allPoolsInOneBSA[[y]]=summedPoolsBSA[[1]][,1:4] #The first four columns are just the information about each mutation, such as CHR and Position...
  cond=grep(y, names(summedPoolsBSA), value=T)
  allPoolsInOneBSA[[y]]$RealDepth=apply(sapply(summedPoolsBSA[cond],function(x) x$RealDepth), 1, sum) #This is another way of doing it (see line 53)
  allPoolsInOneBSA[[y]]$Reference=apply(sapply(summedPoolsBSA[cond],function(x) x$Reference), 1, sum)
  allPoolsInOneBSA[[y]]$Alternative1=apply(sapply(summedPoolsBSA[cond],function(x) x$Alternative1), 1, sum)
}

#This will analyze the comparisons between the conditions of the allPoolsInOne object

for (i in c("Y","L", "B")){
  lowbulk="I"
  highbulk=i
  allPoolsInOneComparisonBSA[[i]]=picker2(allPoolsInOneBSA, lowBulk = lowbulk, highBulk = highbulk)
  print(paste("Comparisons made:",highbulk,"x", lowbulk, sep = " "))
}

###FILTERING

#replicates_filteredBSA=list()
pools_filteredBSA=list()
allPoolsInOne_filteredBSA=list()


for (i in c(names(poolsBSA)[c(-3,-6,-9,-12,-15)],names(poolsBSA)[c(3,6,9,12,15)])){
  print(i)
  pools_filteredBSA[[i]]= analyzer(poolsBSA[[i]][!is.nan(poolsBSA[[i]]$deltaSNP),], minDepthPercentile = 0.005, 
                                   maxDepthPercentile = 0.9, outlierFilt = "deltaSNP")
}

for (i in names(allPoolsInOneComparisonBSA)){
  print(i)
  allPoolsInOne_filteredBSA[[i]]= analyzer(allPoolsInOneComparisonBSA[[i]][!is.nan(allPoolsInOneComparisonBSA[[i]]$deltaSNP),], 
                                           minDepthPercentile = 0.005, maxDepthPercentile = 0.995, outlierFilt = "deltaSNP")
}

save.image(file = "BSA2.RData")
