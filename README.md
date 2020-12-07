# Indicators of inventory completeness for global biodiversity databases
Juliana Stropp<br/>

This document presents the R code used to calculate indicators of inventory completeness presented in Stropp et al. (in prep.).<br/>

R functions were adapted from: https://github.com/AndreMenegotto/SpatialGaps

### Sample coverage (Chao and Jost, 2012) using number of unique sampling events (as in Menegotto & Rangel 2018)
```
ChaoJost_SampEvent <- function(matriz)
    {
  # Determine species unique sampling events
  unic <- unique(matriz[,c(4,5,6,7,8)])
  
  # Define groups to calculate completeness: cell_id, synth_com, sampling_intensity
  #unic[,5] <- cut(unic[,3], breaks = seq(-80, 85, 5))
  #lev <- levels(unic[,5])
  #unic[,5]<-matriz[,c(1,2,3)]
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
### Inventory completeness based on Chao1 using sampling events

```
Chao1_SampEvent <- function(matriz)
{
  # Determine species unique sampling events
  unic <- unique(matriz[,c(4,5,6,7,8)])
  
  # Define groups to calculate completeness: cell_id, synth_com, sampling_intensity
  #unic[,5] <- cut(unic[,3], breaks = seq(-80, 85, 5))
  #lev <- levels(unic[,5])
  #unic[,5]<-matriz[,c(1,2,3)]
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
      # Completeness estimate
      # Number of observed species
      sp <- unique(unic[rec,1])
      Sobs <- length(sp)
      n <- length(rec)
      # Number of observed species only once
      a <- sum(table(unic[rec,1])==1)
      # Number of observed twice
      b <- sum(table(unic[rec,1])==2)
      
      if(b==0)
      {
        k[i] <- NA
      }
      else
      {
        k[i] <- Sobs/(Sobs+((a^2)/(2*b))*((n-1)/n))
      }
    }
  }
  return(k)
}
```
