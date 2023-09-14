
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(parallel)
library(Rmisc)
library(grid)
library(gridExtra)

load(file = "BSA2.RData")

###PLOTS

chosenSize=5000
Pools=list()
allPools=list()

ptm=proc.time()
no_cores <- length(names(pools_filteredBSA))#detectCores()-1 #Use all cores but one, so the computer can do other things while it runs
cl <- makeCluster(no_cores)
clusterExport(cl,c("pools_filteredBSA","allPoolsInOne_filteredBSA" ,"binner2", "Pools","allPools", "chosenSize"))
#Calling the binner2 in parallel
Pools=parLapply(cl,names(pools_filteredBSA),
                function(x){binner2(pools_filteredBSA[[grep(x, names(pools_filteredBSA))]], bin.size = chosenSize)})
#allPools=parLapply(cl,names(allPoolsInOne_filteredBSA),
#                  function(x){binner2(allPoolsInOne_filteredBSA[[grep(x, names(allPoolsInOne_filteredBSA))]], bin.size = chosenSize)})
stopCluster(cl) #Closing the parallel cluster
proc.time()-ptm

for(i in names(allPoolsInOne_filteredBSA)){
  allPools[[i]]=binner2(allPoolsInOne_filteredBSA[[i]], 
                        bin.size = chosenSize)
}

names(Pools)=names(pools_filteredBSA)
names(allPools)=names(allPoolsInOne_filteredBSA)
condition=c("YPD","Lungs")

######For plots with smoothed Delta SNP #####
samples=list()
samples=lapply(c(1:3), FUN = function(x){
  data.frame(CHROM=allPools[[names(allPools)[x]]]$CHROM, 
             binFloor=allPools[[names(allPools)[x]]]$binFloor,
             binCeiling=allPools[[names(allPools)[x]]]$binCeiling,
             binMiddle=allPools[[names(allPools)[x]]]$binMiddle,
             Bin_size=allPools[[names(allPools)[x]]]$Bin_size,
             P1=Pools[[paste0("P1",names(allPools)[x])]]$smoothedDeltaSNP,
             P1sig=Pools[[paste0("P1",names(allPools)[x])]]$Significance,
             P2=Pools[[paste0("P2",names(allPools)[x])]]$smoothedDeltaSNP,
             P2sig=Pools[[paste0("P2",names(allPools)[x])]]$Significance,
             P3=Pools[[paste0("P3",names(allPools)[x])]]$smoothedDeltaSNP,
             P3sig=Pools[[paste0("P3",names(allPools)[x])]]$Significance,
             P4=Pools[[paste0("P4",names(allPools)[x])]]$smoothedDeltaSNP,
             P4sig=Pools[[paste0("P4",names(allPools)[x])]]$Significance,
             P5=Pools[[paste0("P5",names(allPools)[x])]]$smoothedDeltaSNP,
             P5sig=Pools[[paste0("P5",names(allPools)[x])]]$Significance,
             All_pools=allPools[[names(allPools)[x]]]$smoothedDeltaSNP,
             mut_perBin=allPools[[names(allPools)[x]]]$mut_perBin,
             Significance=allPools[[names(allPools)[x]]]$Significance) %>%
    mutate(pool_max=pmax(P1,P2,P3,P4,P5, na.rm=T), 
           pool_min=pmin(P1,P2,P3,P4,P5, na.rm=T), 
           Condition=rep(condition[x], nrow(allPools[[names(allPools)[x]]])) ) 
})
names(samples)=names(allPools)[c(1:3)]

######For plots with NON-smoothed Delta SNP #####
samplesNonSmooth=list()
samplesNonSmooth=lapply(c(1:3), FUN = function(x){
  data.frame(CHROM=allPools[[names(allPools)[x]]]$CHROM, 
             binFloor=allPools[[names(allPools)[x]]]$binFloor,
             binCeiling=allPools[[names(allPools)[x]]]$binCeiling,
             binMiddle=allPools[[names(allPools)[x]]]$binMiddle,
             Bin_size=allPools[[names(allPools)[x]]]$Bin_size,
             P1=Pools[[paste0("P1",names(allPools)[x])]]$smoothedDeltaSNP,
             P1sig=Pools[[paste0("P1",names(allPools)[x])]]$Significance,
             P2=Pools[[paste0("P2",names(allPools)[x])]]$smoothedDeltaSNP,
             P2sig=Pools[[paste0("P2",names(allPools)[x])]]$Significance,
             P3=Pools[[paste0("P3",names(allPools)[x])]]$smoothedDeltaSNP,
             P3sig=Pools[[paste0("P3",names(allPools)[x])]]$Significance,
             P4=Pools[[paste0("P4",names(allPools)[x])]]$smoothedDeltaSNP,
             P4sig=Pools[[paste0("P4",names(allPools)[x])]]$Significance,
             P5=Pools[[paste0("P5",names(allPools)[x])]]$smoothedDeltaSNP,
             P5sig=Pools[[paste0("P5",names(allPools)[x])]]$Significance,

             All_pools=allPools[[names(allPools)[x]]]$deltaSNP,
             mut_perBin=allPools[[names(allPools)[x]]]$mut_perBin,
             Significance=allPools[[names(allPools)[x]]]$Significance) %>%
    mutate(pool_max=pmax(P1,P2,P3,P4,P5, na.rm=T), 
           pool_min=pmin(P1,P2,P3,P4,P5, na.rm=T), 
           Condition=rep(condition[x], nrow(allPools[[names(allPools)[x]]])) ) 
})

names(samplesNonSmooth)=names(allPools)[c(1:3)]




#This is a table with the Max and Min of ∆SNP for each position.
MaxAndMins=cbind(allPoolsInOne_filteredBSA[["Y"]][,1:4], data.frame(Min_SNPchange=-(allPoolsInOne_filteredBSA[["Y"]]$SNPindex.LOW),
                                                                    Max_SNPchange=1-(allPoolsInOne_filteredBSA[["Y"]]$SNPindex.LOW))  ) #The SNPindex.LOW should be the same for YPD and Lungs.... 


MaxAndMins=MaxAndMins %>% filter(!(CHROM=="chr2" & (POS==283162 | POS==283652| (POS>466000 & POS<467000)   ) )) #This mutation was erased in the addition of the drug marker in the C8 parent (It was in the arm).

#PLOTS
#The IRs
IR_list=list()
IR_list[["chr1"]]=c(1262.121, 1340.862)
IR_list[["chr2"]]=c(283.162,	468.855)
IR_list[["chr3"]]=c(15.943 , 117.260)
IR_list[["chr7"]]=c(19.297, 146.051)
IR_list[["chr8"]]=c(13.998, 55.13)


# Defining temp, a subset with only the chosen chromosome.

dataToUse=samplesNonSmooth

for(i in 1:14){
  rm(temp)
  rm(PoolsCombine)
  rm(newMinsAndMax)
  
  temp=rbind(dataToUse[["L"]][dataToUse[["L"]]$CHROM==paste0("chr",i),], 
             dataToUse[["Y"]][dataToUse[["Y"]]$CHROM==paste0("chr",i),])
  temp$combine=paste(temp$Significance, temp$Condition, sep="_")
  levels(temp$combine)=c("FALSE_Inoculum", "FALSE_YPD", "TRUE_Inoculum","TRUE_YPD")
  eachPool=c("P1sig","P2sig","P3sig","P4sig","P5sig")
  PoolsCombine=lapply(eachPool, FUN=function(x){z=paste(temp[,grep(x,colnames(temp))], temp$Condition, sep="_");return(z)   })
  
  newMinsAndMax= MaxAndMins %>% filter(CHROM== paste0("chr",i))
  
  p1=ggplot(data=temp,aes(x=temp$binMiddle/1000, y=temp$All_pools))+
    geom_line(data = newMinsAndMax, aes(x=POS/1000, y=Min_SNPchange))+
    geom_line(data = newMinsAndMax, aes(x=POS/1000, y=Max_SNPchange))+
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
    scale_y_continuous(breaks = c(-1,-0.75,-0.5, -0.25, 0,0.25,0.5,0.75,1), limits = c(-1.05,1.05), expand = c(0,0))+
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
  
   ggsave(p1,filename = paste0("/Users/danielagustinho/GoogleDrive/WashU/DNA/EXP#067/Plots/Chromosome_", i, "_nonSmoothed.1kb.pdf"), dpi=500,height = 9,
         units = "in", width = 16)
  
  ggsave(p2,filename = paste0("/Users/danielagustinho/GoogleDrive/WashU/DNA/EXP#067/Plots/Chromosome_", i, "_Lungs15days_nonSmoothed.1kb.pdf"), dpi=500,height = 9,
         units = "in", width = 16)
  
  ggsave(p3,filename = paste0("/Users/danielagustinho/GoogleDrive/WashU/DNA/EXP#067/Plots/Chromosome_", i, "_Brains15days_nonSmoothed.1kb.pdf"), dpi=500,height = 9,
         units = "in", width = 16)
}


####IR-1 PLOTS
i=2
rm(temp)
rm(PoolsCombine)
rm(newMinsAndMax)
temp=rbind(samplesNonSmooth[["L"]][samplesNonSmooth[["L"]]$CHROM==paste0("chr",i),], 
           samplesNonSmooth[["Y"]][samplesNonSmooth[["Y"]]$CHROM==paste0("chr",i),])
temp$combine=paste(temp$Significance, temp$Condition, sep="_")
levels(temp$combine)=c("FALSE_Lungs", "FALSE_YPD", "TRUE_Lungs","TRUE_YPD")
eachPool=c("P1sig","P2sig","P3sig","P4sig","P5sig")
PoolsCombine=lapply(eachPool, FUN=function(x){z=paste(temp[,grep(x,colnames(temp))], temp$Condition, sep="_");return(z)   })

newMinsAndMax= MaxAndMins %>% filter(CHROM== paste0("chr",i))

p1=ggplot(data=temp,aes(x=temp$binMiddle/1000, y=temp$All_pools))+
  geom_rect(xmin  =283.162, xmax=	468.855,ymin=-1.1, ymax=1.1,color="yellow",fill="yellow",alpha=0.01)+
  geom_line(data = newMinsAndMax, aes(x=POS/1000, y=Min_SNPchange))+
  geom_line(data = newMinsAndMax, aes(x=POS/1000, y=Max_SNPchange))+
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
  labs(title=paste0("Chromosome ", i), x=paste0("Position on chromosome ",i," (kb)"),y="Change in allele frequency (BSA)", color="Pool")+
  geom_hline(yintercept = c(0), color="black", size=0.6)+
  scale_x_continuous(breaks = seq(from=0, to=2400, by=300), limits=c(0, 1650),expand = c(0, 0))+
  scale_y_continuous(breaks = c(-1,-0.75, -0.5, -0.25, 0,0.25,0.5,0.75,1), limits = c(-1.1,1.1),expand = c(0, 0) )+
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
p1

ggsave(p1,filename = paste0("/Users/danielagustinho/GoogleDrive/WashU/DNA/EXP#067/Plots/IR-1_5kb-NONsmoothed.pdf"), dpi=500,height = 9,
       units = "in", width = 16)

###IR-1 Zoom PLOT
PoolsZoom=list()
allPoolsZoom=list()
for(i in names(pools_filteredBSA)){
  PoolsZoom[[i]]=pools_filteredBSA[[i]] %>% dplyr::filter( (CHROM=="chr2" & POS>249500  & POS<500500))%>%
    mutate(Significance=F)
  for(j in 1:nrow(PoolsZoom[[i]])){if(PoolsZoom[[i]]$qvalue[j]<0.1){PoolsZoom[[i]]$Significance[j]=TRUE}}
}
for(i in names(allPoolsInOne_filteredBSA)){
  allPoolsZoom[[i]]=allPoolsInOne_filteredBSA[[i]]%>% dplyr::filter( (CHROM=="chr2" & POS>249500  & POS<500500))%>%
    mutate(Significance=F)
  for(j in 1:nrow(allPoolsZoom[[i]])){if(allPoolsZoom[[i]]$qvalue[j]<0.1){allPoolsZoom[[i]]$Significance[j]=TRUE}}
}


names(PoolsZoom)=c("P1.Y", "P1.L","P2.Y","P2.L","P3.Y","P3.L","P4.Y","P4.L","P5.Y","P5.L","Pall.Y","Pall.L")
names(allPoolsZoom)=c("Y", "L")
condition=c("YPD","Inoculum")

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
condition=c("YPD","Inoculum")

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
levels(temp$combine)=c("FALSE_Inoculum", "FALSE_YPD", "TRUE_Inoculum","TRUE_YPD")
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
                      values = c('FALSE_Inoculum'="red",'FALSE_YPD'="black",
                                 'TRUE_Inoculum'="red",'TRUE_YPD'="black"), drop=F)+
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


ggsave(p1,filename = paste0("/Users/danielagustinho/GoogleDrive/WashU/DNA/EXP#067/Plots/IR-1zoom.pdf"), dpi=500,height = 9,units = "in", width = 16)

