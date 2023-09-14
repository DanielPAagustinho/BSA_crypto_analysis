
library(QTLseqr)
library(ggplot2)
library(stringr)
library(tidyr)
library(reshape2)
library(RColorBrewer)
library(parallel)
library(Rmisc)
library(grid)
library(gridExtra)
library(tidyverse)

source("BSA6tools.R")

#####Parent strains -
parents=list()
for(x in c("KN99a", "TDY1993")){
  parents[[x]]=read.table(file=paste0("Parents_Tables/",x,".txt"),header = T,stringsAsFactors = F, numerals = "allow.loss");
  parents[[x]]$Alternative1=as.integer(parents[[x]]$Alternative1)
}
#Finding the positions where KN99a is reference and TDY1993 (C8) is alternative.
parents=parents[[1]][parents[[2]]$Filt_Genotype=="Alternative" & parents[[1]]$Filt_Genotype== "Reference",] 





poolsBSA4=list.files(path="VCF_tables") %>%
  sub(pattern=".txt", replacement = "")
mySamples=list()
for (i in poolsBSA4){
  mySamples[[i]]=read.table(file = paste0("VCF_tables/",i,
                                          ".txt"), header = T,stringsAsFactors = F, numerals = "allow.loss");
  mySamples[[i]]=mySamples[[i]][paste(mySamples[[i]]$CHR, mySamples[[i]]$POS, sep="_") %in% paste(parents$CHR, parents$POS, sep="_") ,]
  mySamples[[i]]$Alternative1=as.integer(mySamples[[i]]$Alternative1)
}

#Summary of data from the coverage
SampleMean=sapply(mySamples, FUN=function(x){mean(x$RealDepth)})
SampleMedian=sapply(mySamples, FUN=function(x){median(x$RealDepth)})

Inoculum=grep(pattern="Pool_[[:alnum:]]{1,3}.C", poolsBSA4, value=T)
YPD=grep(pattern="Pool_[[:alnum:]]{1,3}.Y", poolsBSA4, value=T)
Lungs=grep(pattern="Pool_[[:alnum:]]{1,3}.[[:digit:]].L", poolsBSA4, value=T)

replicatesBSA4=list()


#Filling up the replicates list.
conditions=c("Y","L") #Conditions to be used as high bulk

for(i in c(1:5, "ALL")){ #The pools
  for(z in 1:3){ # The replicates
    for(y in 1:length(conditions)){
      lowbulk=paste0("Pool_", i, ".","C")
      if (y=="1"){highbulk=paste0("Pool_", i, ".",conditions[y])} else {highbulk=paste0("Pool_", i, ".",z,".",conditions[y])}
      if(lowbulk %in% names(mySamples) & highbulk %in% names(mySamples) ){
        replicatesBSA4[[highbulk]]=picker2(mySamples, lowBulk = lowbulk, highBulk = highbulk)
        print(paste("Comparisons made:",lowbulk,"x", highbulk, sep = " "))}
    }
  }
}


#Making a new list with the sums of all replicates of a pool

summedPoolsBSA4=list() # a list with the summed depth and allele counts for all replicates of a specific pool
poolsBSA4=list() # a list of comparisons of pooled replicates (all the replicates from a specific pool) for each condition using picker()
listOfPoolsBSA4=gsub(".C","",Inoculum)
ConditionsBSA4=c("L") # The one condition that was done in triplicate: Lungs

for(i in listOfPoolsBSA4){ #This loop create the table that is the sum of all replicates of a pool.
  for(y in ConditionsBSA4){
    reps=grep(paste0(i, "...",y),names(mySamples))
    meta=mySamples[[reps[1]]][,1:4]
    if(length(grep(paste0(i, "...",y),names(mySamples)))==3){
      meta$RealDepth=mySamples[[reps[1]]]$RealDepth+mySamples[[reps[2]]]$RealDepth+mySamples[[reps[3]]]$RealDepth
      meta$Reference=mySamples[[reps[1]]]$Reference+mySamples[[reps[2]]]$Reference+mySamples[[reps[3]]]$Reference
      meta$Alternative1=mySamples[[reps[1]]]$Alternative1+mySamples[[reps[2]]]$Alternative1+mySamples[[reps[3]]]$Alternative1
    }
    if(length(grep(paste0(i, "...",y),names(mySamples)))==2){
      meta$RealDepth=mySamples[[reps[1]]]$RealDepth+mySamples[[reps[2]]]$RealDepth
      meta$Reference=mySamples[[reps[1]]]$Reference+mySamples[[reps[2]]]$Reference
      meta$Alternative1=mySamples[[reps[1]]]$Alternative1+mySamples[[reps[2]]]$Alternative1
    }
    summedPoolsBSA4[[paste0(i,".",y)]]=meta
    rm(meta)
  }
 }
for(j in c(Inoculum,YPD)){summedPoolsBSA4[[j]]=mySamples[[j]]}

##Making the contrasts of low and high bulk with the summed pools

for(i in listOfPoolsBSA4){ #The pools
  for(y in conditions){
    lowbulk=paste0(i,".C")
    highbulk=paste0(i, ".",y)
    if(lowbulk %in% names(summedPoolsBSA4) & highbulk %in% names(summedPoolsBSA4) ){
      poolsBSA4[[highbulk]]=picker2(summedPoolsBSA4, lowBulk = lowbulk, highBulk = highbulk)
      print(paste("Comparisons made:",lowbulk,"x", highbulk, sep = " "))}
    
  }
}

#Here we start summing all pools together for each condition, and then making the comparisons.

allPoolsInOneBSA4=list() #This will be the list containing three itens, the sum of all inocula, all YPDs and all lungs.
allPoolsInOneComparisonBSA4=list()

for (y in c("C","Y",ConditionsBSA4)){ #This loop goes around for each condition.
  allPoolsInOneBSA4[[y]]=summedPoolsBSA4[[1]][,1:4] #The first four columns are just the information about each mutation, such as CHR and Position...
  cond=grep(paste0("Pool_","[[:alnum:]]{1,3}.",y), names(summedPoolsBSA4), value=T)[1:6] #1 to 6 because there are only 6 different pools...
  allPoolsInOneBSA4[[y]]$RealDepth=apply(sapply(summedPoolsBSA4[cond],function(x) x$RealDepth), 1, sum)
  allPoolsInOneBSA4[[y]]$Reference=apply(sapply(summedPoolsBSA4[cond],function(x) x$Reference), 1, sum)
  allPoolsInOneBSA4[[y]]$Alternative1=apply(sapply(summedPoolsBSA4[cond],function(x) x$Alternative1), 1, sum)
}

#This will analyze the comparisons between the conditions of the allPoolsInOne object

for (i in conditions){
  lowbulk="C"
  highbulk=i
  allPoolsInOneComparisonBSA4[[i]]=picker2(allPoolsInOneBSA4, lowBulk = lowbulk, highBulk = highbulk)
  print(paste("Comparisons made:",lowbulk,"x", highbulk, sep = " "))
}

###FILTERING

#replicates_filteredBSA4=list()
pools_filteredBSA4=list()
allPoolsInOne_filteredBSA4=list()


for (i in names(poolsBSA4)){
  print(i)
  pools_filteredBSA4[[i]]= analyzer(poolsBSA4[[i]][!is.nan(poolsBSA4[[i]]$deltaSNP),], minDepthPercentile = 0.1, maxDepthPercentile = 0.995)
}

for (i in names(allPoolsInOneComparisonBSA4)){
  print(i)
  allPoolsInOne_filteredBSA4[[i]]= analyzer(allPoolsInOneComparisonBSA4[[i]][!is.nan(allPoolsInOneComparisonBSA4[[i]]$deltaSNP),], 
                                            minDepthPercentile = 0.1, maxDepthPercentile = 0.995)
}


###PLOTS

chosenSize=10000
Pools=list()
allPools=list()

ptm=proc.time()
no_cores <- detectCores()-1 #Use all cores but one, so the computer can do other things while it runs
cl <- makeCluster(no_cores)
clusterExport(cl,c("pools_filteredBSA4","allPoolsInOne_filteredBSA4" ,"binner2", "Pools","allPools", "chosenSize"))
#Calling the binner2 in parallel
Pools=parLapply(cl,names(pools_filteredBSA4),
                function(x){binner2(pools_filteredBSA4[[grep(x, names(pools_filteredBSA4))]], bin.size = chosenSize)})
allPools=parLapply(cl,names(allPoolsInOne_filteredBSA4),
                   function(x){binner2(allPoolsInOne_filteredBSA4[[grep(x, names(allPoolsInOne_filteredBSA4))]], bin.size = chosenSize)})
stopCluster(cl) #Closing the parallel cluster
proc.time()-ptm

names(Pools)=c("P1.Y", "P1.L","P2.Y","P2.L","P3.Y","P3.L","P4.Y","P4.L","P5.Y","P5.L","Pall.Y","Pall.L")
names(allPools)=c("Y", "L")
condition=c("YPD","Lungs")

samples=lapply(c(1,2), FUN = function(x){
  data.frame(CHROM=allPools[[names(allPools)[x]]]$CHROM, 
             binFloor=allPools[[names(allPools)[x]]]$binFloor,
             binCeiling=allPools[[names(allPools)[x]]]$binCeiling,
             binMiddle=allPools[[names(allPools)[x]]]$binMiddle,
             Bin_size=allPools[[names(allPools)[x]]]$Bin_size,
             P1=Pools[[paste0("P1.",names(allPools)[x])]]$smoothedDeltaSNP,
             P1sig=Pools[[paste0("P1.",names(allPools)[x])]]$Significance,
             P2=Pools[[paste0("P2.",names(allPools)[x])]]$smoothedDeltaSNP,
             P2sig=Pools[[paste0("P2.",names(allPools)[x])]]$Significance,
             P3=Pools[[paste0("P3.",names(allPools)[x])]]$smoothedDeltaSNP,
             P3sig=Pools[[paste0("P3.",names(allPools)[x])]]$Significance,
             P4=Pools[[paste0("P4.",names(allPools)[x])]]$smoothedDeltaSNP,
             P4sig=Pools[[paste0("P4.",names(allPools)[x])]]$Significance,
             P5=Pools[[paste0("P5.",names(allPools)[x])]]$smoothedDeltaSNP,
             P5sig=Pools[[paste0("P5.",names(allPools)[x])]]$Significance,
             Pall=Pools[[paste0("Pall.",names(allPools)[x])]]$smoothedDeltaSNP,
             Pallsig=Pools[[paste0("Pall.",names(allPools)[x])]]$Significance,
             All_pools=allPools[[names(allPools)[x]]]$smoothedDeltaSNP,
             mut_perBin=allPools[[names(allPools)[x]]]$mut_perBin,
             Significance=allPools[[names(allPools)[x]]]$Significance) %>%
    mutate(pool_max=pmax(P1,P2,P3,P4,P5, na.rm=T), pool_min=pmin(P1,P2,P3,P4,P5, na.rm=T), Condition=rep(condition[x], nrow(allPools[[names(allPools)[x]]]))   ) 
})

names(samples)=names(allPools)


save.image(file = "GoogleDrive/WashU/DNA/EXP#054/Rdata.RData")
load(file = "GoogleDrive/WashU/DNA/EXP#054/Rdata.RData")


#PLOTS
#The IRs
IR_list=list()
IR_list[["chr1"]]=c(1262.121, 1340.862)
IR_list[["chr2"]]=c(283.162,	468.855)
IR_list[["chr3"]]=c(15.943 , 117.260)
IR_list[["chr7"]]=c(19.297, 146.051)
IR_list[["chr8"]]=c(13.998, 55.13)


# Defining temp, a subset with only the chosen chromosome.
for(i in 1:14){
  rm(temp)
  rm(PoolsCombine)
  temp=rbind(samples[["L"]][samples[["L"]]$CHROM==paste0("chr",i),], 
             samples[["Y"]][samples[["Y"]]$CHROM==paste0("chr",i),])
  temp$combine=paste(temp$Significance, temp$Condition, sep="_")
  levels(temp$combine)=c("FALSE_Lungs", "FALSE_YPD", "TRUE_Lungs","TRUE_YPD")
  eachPool=c("P1sig","P2sig","P3sig","P4sig","P5sig", "Pallsig")
  PoolsCombine=lapply(eachPool, FUN=function(x){z=paste(temp[,grep(x,colnames(temp))], temp$Condition, sep="_");return(z)   })

p1=ggplot(data=temp,aes(x=temp$binMiddle/1000, y=temp$All_pools))+
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


ggsave(p1,filename = paste0("GoogleDrive/WashU/DNA/EXP#054/Plots_grant/allChromosomes_chr", i, ".10kb.pdf"), dpi=500,height = 9,
       units = "in", width = 16)

}


####IR-1 PLOTS
i=2
rm(temp)
rm(PoolsCombine)
temp=rbind(samples[["L"]][samples[["L"]]$CHROM==paste0("chr",i),], 
           samples[["Y"]][samples[["Y"]]$CHROM==paste0("chr",i),])
temp$combine=paste(temp$Significance, temp$Condition, sep="_")
levels(temp$combine)=c("FALSE_Lungs", "FALSE_YPD", "TRUE_Lungs","TRUE_YPD")
eachPool=c("P1sig","P2sig","P3sig","P4sig","P5sig", "Pallsig")
PoolsCombine=lapply(eachPool, FUN=function(x){z=paste(temp[,grep(x,colnames(temp))], temp$Condition, sep="_");return(z)   })

p1=ggplot(data=temp,aes(x=temp$binMiddle/1000, y=temp$All_pools))+
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
  geom_vline(xintercept  = c(283.162,	468.855),color="black", size=0.3, linetype="dashed")+
  #geom_vline(xintercept = 0, size=0.3)+
  labs(title=paste0("Chromosome ", i), x=paste0("Position on chromosome ",i," (kb)"),y="Change in allele frequency (BSA)", color="Pool")+
  geom_hline(yintercept = c(0), color="black", size=0.6)+
  scale_x_continuous(breaks = seq(from=0, to=2400, by=300), limits=c(0, 1650),expand = c(0, 0))+
  scale_y_continuous(breaks = c(-0.5, -0.25, 0,0.25,0.5), limits = c(-0.60,0.60), expand = c(0,0))+
  theme(#legend.title = element_text(size=18), 
    legend.position = "none", 
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


ggsave(p1,filename = paste0("GoogleDrive/WashU/DNA/EXP#054/Plots_grant/IR-1_10kb.pdf"), dpi=500,height = 9,
       units = "in", width = 16)

###IR-1 Zoom PLOT
PoolsZoom=list()
allPoolsZoom=list()
for(i in names(pools_filteredBSA4)){
  PoolsZoom[[i]]=pools_filteredBSA4[[i]] %>% dplyr::filter( (CHROM=="chr2" & POS>249500  & POS<500500))%>%
    mutate(Significance=F)
  for(j in 1:nrow(PoolsZoom[[i]])){if(PoolsZoom[[i]]$qvalue[j]<0.1){PoolsZoom[[i]]$Significance[j]=TRUE}}
}
for(i in names(allPoolsInOne_filteredBSA4)){
  allPoolsZoom[[i]]=allPoolsInOne_filteredBSA4[[i]]%>% dplyr::filter( (CHROM=="chr2" & POS>249500  & POS<500500))%>%
    mutate(Significance=F)
  for(j in 1:nrow(allPoolsZoom[[i]])){if(allPoolsZoom[[i]]$qvalue[j]<0.1){allPoolsZoom[[i]]$Significance[j]=TRUE}}
}


names(PoolsZoom)=c("P1.Y", "P1.L","P2.Y","P2.L","P3.Y","P3.L","P4.Y","P4.L","P5.Y","P5.L","Pall.Y","Pall.L")
names(allPoolsZoom)=c("Y", "L")
condition=c("YPD","Lungs")

chosenSize=500
Pools=list()
allPools=list()

ptm=proc.time()
no_cores <- detectCores()-1 #Use all cores but one, so the computer can do other things while it runs
cl <- makeCluster(no_cores)
clusterExport(cl,c("PoolsZoom","allPoolsZoom" ,"binner2", "Pools","allPools", "chosenSize"))
#Calling the binner2 in parallel
Pools=parLapply(cl,names(PoolsZoom),
                function(x){binner2(PoolsZoom[[grep(x, names(PoolsZoom))]], bin.size = chosenSize)})
allPools=parLapply(cl,names(allPoolsZoom),
                   function(x){binner2(allPoolsZoom[[grep(x, names(allPoolsZoom))]], bin.size = chosenSize)})
stopCluster(cl) #Closing the parallel cluster
proc.time()-ptm

names(Pools)=c("P1.Y", "P1.L","P2.Y","P2.L","P3.Y","P3.L","P4.Y","P4.L","P5.Y","P5.L","Pall.Y","Pall.L")
names(allPools)=c("Y", "L")
condition=c("YPD","Lungs")

samples=lapply(c(1,2), FUN = function(x){
  data.frame(CHROM=allPools[[names(allPools)[x]]]$CHROM, 
             binFloor=allPools[[names(allPools)[x]]]$binFloor,
             binCeiling=allPools[[names(allPools)[x]]]$binCeiling,
             binMiddle=allPools[[names(allPools)[x]]]$binMiddle,
             Bin_size=allPools[[names(allPools)[x]]]$Bin_size,
             P1=Pools[[paste0("P1.",names(allPools)[x])]]$smoothedDeltaSNP,
             P1sig=Pools[[paste0("P1.",names(allPools)[x])]]$Significance,
             P2=Pools[[paste0("P2.",names(allPools)[x])]]$smoothedDeltaSNP,
             P2sig=Pools[[paste0("P2.",names(allPools)[x])]]$Significance,
             P3=Pools[[paste0("P3.",names(allPools)[x])]]$smoothedDeltaSNP,
             P3sig=Pools[[paste0("P3.",names(allPools)[x])]]$Significance,
             P4=Pools[[paste0("P4.",names(allPools)[x])]]$smoothedDeltaSNP,
             P4sig=Pools[[paste0("P4.",names(allPools)[x])]]$Significance,
             P5=Pools[[paste0("P5.",names(allPools)[x])]]$smoothedDeltaSNP,
             P5sig=Pools[[paste0("P5.",names(allPools)[x])]]$Significance,
             Pall=Pools[[paste0("Pall.",names(allPools)[x])]]$smoothedDeltaSNP,
             Pallsig=Pools[[paste0("Pall.",names(allPools)[x])]]$Significance,
             All_pools=allPools[[names(allPools)[x]]]$smoothedDeltaSNP,
             mut_perBin=allPools[[names(allPools)[x]]]$mut_perBin,
             Significance=allPools[[names(allPools)[x]]]$Significance) %>%
    mutate(pool_max=pmax(P1,P2,P3,P4,P5, na.rm=T), pool_min=pmin(P1,P2,P3,P4,P5, na.rm=T), Condition=rep(condition[x], nrow(allPools[[names(allPools)[x]]]))   ) 
})

names(samples)=names(allPools)
samplesZoom=samples

i=2
rm(temp)
rm(PoolsCombine)
temp=rbind(samplesZoom[["L"]][samplesZoom[["L"]]$CHROM==paste0("chr",i),], 
           samplesZoom[["Y"]][samplesZoom[["Y"]]$CHROM==paste0("chr",i),])
temp$combine=paste(temp$Significance, temp$Condition, sep="_")
levels(temp$combine)=c("FALSE_Lungs", "FALSE_YPD", "TRUE_Lungs","TRUE_YPD")
eachPool=c("P1sig","P2sig","P3sig","P4sig","P5sig", "Pallsig")
PoolsCombine=lapply(eachPool, FUN=function(x){z=paste(temp[,grep(x,colnames(temp))], temp$Condition, sep="_");return(z)   })

p1=ggplot(data=temp,aes(x=temp$binMiddle/1000, y=temp$All_pools))+
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
  geom_vline(xintercept  = c(283.162,	468.855),color="black", size=0.3, linetype="dashed")+
  #geom_vline(xintercept = 0, size=0.3)+
  labs(title=paste0("Chromosome ", i), x=paste0("Position on chromosome ",i," (kb)"),y="Change in allele frequency (BSA)", color="Pool")+
  geom_hline(yintercept = c(0), color="black", size=0.6)+
  scale_x_continuous(breaks = seq(from=0, to=2400, by=50), limits=c(250, 500),expand = c(0, 0))+
  scale_y_continuous(breaks = c(-0.5, -0.25, 0,0.25,0.5), limits = c(-0.60,0.60), expand = c(0,0))+
  theme(#legend.title = element_text(size=18), 
    legend.position = "none", 
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


ggsave(p1,filename = paste0("GoogleDrive/WashU/DNA/EXP#054/Plots_grant/IR-1zoom.pdf"), dpi=500,height = 9,units = "in", width = 16)
       
       
       
      