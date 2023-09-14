# The tabler() funtion that I used above is defined below. You'll need to install the "stringr" package on R to make it work: 
# install.packages("stringr").

#  This function extracts info from an imported vcf for 1 sample, outputting a nice table.

#  @param VCF 
#  @param i integer. Samples in a vcf file start in the 10th column. Tabler takes the information from the (9+i)th column from the vcf.
#  @param minDepth Numeric. The minimal depth required at each position to call a different allele. Default=5.
#  @param cutoff Numeric. Minimum ratio of reads (as a number between 0 and 1) that supports either allele required to call an allele. Default=0.75. 
#  @param filtering Logical. If True, will fill up the column "filtGenotype". This is specially useful for single strains, not so much for BSA. The function is much faster if turned off. Default=F. 
#  @return a table

tabler = function(VCF, i=1, minDepth=5, cutoff=0.75, filtering=F)
{
library(stringr)
  c=(str_split_fixed(VCF[,9+i],":",n=8)) #This splits the FORMAT items from the vcf into columns.
  colnames(c)=c("Genotype_GT","Depth_DP","AD","Reference_RO","QR","Alternative_AO","QA","GL")
  #The two lines below convert any "." in the object c (columns 3, 8 and the rest, respectively) to "0,0".
  c[,3]=sapply(c[,3], USE.NAMES = F, FUN = function(x){if(x=="."){x="0,0"} else {x=x}  })
  c[,8]=sapply(c[,8], USE.NAMES = F, FUN = function(x){if(x=="."){x="0,0"} else {x=x}  })
  c=apply(c,c(1,2), FUN = function(x){if(x=="."){x=0} else {x=x}  }  )
  
  #This gets the the number of Reference reads and Alternative reads from the column AD from the object c.
  d=str_split_fixed(c[,3], ",",n=3) #This splits the two alternative allele counts into different columns.
  
  #Convert c and d to numeric
  c=apply(c,c(1,2),as.numeric) #Columns AD and GL will be converted to NA, generating errors. But that's ok, we've got what we needed from them. 
  d=apply(d,c(1,2),as.numeric)
  
  #This gets the possible alternative alleles. The split only happens when we have a second alternative allele.
  e=str_split_fixed(VCF[,5], ",",n=3)

  data2=data.frame(CHR=VCF[,1], POS=VCF[,2], Ref_allele= VCF[,4], ALT1=e[,1],ALT2=e[,2],ALT3=e[,3], Genotype=as.numeric(c[,1]),
                   Depth=as.numeric(c[,2]),Reference= as.numeric(c[,4]), Alternative1=d[,2], Alternative2=d[,3], 
                   filtGenotype=rep(".",nrow(e)), stringsAsFactors = FALSE) 
  data2$Alternative2=sapply(data2$Alternative2, USE.NAMES = F, FUN = function(x){if(is.na(x)){x=as.numeric(0)} else {x=x}  })
  
  data2$filtGenotype=as.character(data2$filtGenotype)
  data2=mutate(data2, realDepth=Reference+Alternative1+Alternative2,
         Ref_percentage=Reference/realDepth,
         Alt1_percentage=Alternative1/realDepth,
         Alt2_percentage=Alternative2/realDepth)
  if(filtering==T){
  for(j in 1:nrow(data2)){ ## Loop for defining the filtered Genotype. Filter by minimum depth and Genotype cutoff. The cutoff is the
    # minimum percentage of reads necessary to define a genotype. Default = 0.75.
    if(is.na(data2$Reference[j])){x="noReads"}
    else if(data2[j,8]< (minDepth+1)){data2[j,12]="lowDepth"}
    else if(data2[j,8]> minDepth & data2$Reference[j]>cutoff*data2$realDepth[j]){data2[j,12]="Reference"}
    else if(data2[j,8]>minDepth & data2$Alternative1[j]>cutoff*data2$realDepth[j]){data2[j,12]="Alternative1"}
    else if(data2[j,8]>minDepth & data2$Alternative2[j]>cutoff*data2$realDepth[j]){data2[j,12]="Alternative2"}
    else if(data2[j,8]>minDepth & (data2$Alternative2[j]<=cutoff*data2$realDepth[j] &  data2$Alternative1[j]<=cutoff*data2$realDepth[j] &
                                   data2$Reference[j]<=cutoff*data2$realDepth[j])  ){data2[j,12]="undefinedGenotype"}
    
  		}
	}
  rm(c)
  rm(d)
  return(data2)
} 
