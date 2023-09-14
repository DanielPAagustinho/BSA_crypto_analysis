getStrainNames <- function(rawStrains)
{
  print(rawStrains)
  strainsList <- list()
  i = 1
  for (strain in rawStrains)
  {
    if (startsWith(strain,"sp")[1]){
      strainsList[[i]] <- paste("sp",substr(strain, nchar(strain)-3+1, nchar(strain)),sep="")
    }
    else{
      strainsList[[i]] <- strain
    }
    i = i+1
  }
  return(strainsList)
}

makeGWASData <- function(rawData,strains)
{
  numCols = 16
  numStrains = length(rawData)/numCols
  finalList <- list()
  for (i in 1:numStrains){
    tempDF = rawData[,((i-1)*numCols+1):(i*numCols)]
    rawName <- colnames(tempDF)[1]
    strainName <- stringr::word(rawName,1,sep = stringr::fixed(".CHR"))
    if (strainName %in% strains)
    {
      finalList[[strainName]] <- tempDF
    }
  }
  return(finalList)
}

buildGenotypeMatrix <- function(strains,processedData)
{
  genotypeMatrix <- matrix(, nrow = length(strains), ncol = dim(rawData)[1])
  counter = 1
  for (strain in strains){
    colName = paste(strain,"filtGenotype",sep=".")
    genotype = processedData[[strain]][[colName]]
    genotypeMatrix[counter,] = genotype
    counter = counter + 1
  }
  
  return(genotypeMatrix)
}

buildVariantMatrix <- function(virulence,genotypeData)
{
  variantMatrix <- matrix(, nrow = dim(virulence)[1], ncol = 2)
  variantMatrix[,1] = virulence[,2]
  variantMatrix[,2] = genotypeData
  filtered = variantMatrix[(variantMatrix[,2]=="Reference" | variantMatrix[,2]=="Alternative1"),]
  
  return(filtered)
}