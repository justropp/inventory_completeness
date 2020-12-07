# Indicators of inventory completeness for global biodiversity databases
This document presents the R code used to calculate indicators of inventory completeness presented in Stropp et al. (in prep.).

Calculate sample coverage (Chao and Jost, 2012): number of singletons and doubletons.<br/>R functions were adapted from: https://github.com/AndreMenegotto/SpatialGaps

## Sample coverage
```
SampCovComplet_speciesName <- function(matriz)
    {
  # Determine species unique sampling events
  unic <- unique(matriz[,c(4,5,6,7,8)])
  
  # Define groups to calculate completeness: cell_id, synth_com, sampling_intensity
  lev <-unique(matriz[,8])

  # Calculate level of completeness
  k <- numeric()
  
  # For each group:
  for(i in 1:length(lev))
  {
    # Find all records
    rec <- which(unic[,5]==lev[i])
    
    if(length(rec)==0)
    {
      k[i] <- NA
    }
    else
    {
      # n = number of specimens (number of records)
      n <- length(rec)
      # f1 = number of singletons
      f1 <- sum(table(unic[rec,1])==1)
      # f2 = number of doubletons
      f2 <- sum(table(unic[rec,1])==2)
      
      # Completeness estimate
      num <- (n-1)*f1
      den <- ((n-1)*f1)+(2*f2)
      
      if(den==0)
      {
        k[i] <- NA
      }
      else
      {
        k[i] <- 1 - ((f1/n)*(num/den))
      }
    }
  }
  return(k)
}
```
