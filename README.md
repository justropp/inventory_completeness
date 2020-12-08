# Indicators of inventory completeness for global biodiversity databases
Juliana Stropp<br/>

This document presents the R code used to calculate indicators of inventory completeness presented in Stropp et al. (in prep.; Table 1).<br/>

R functions were adapted from: https://github.com/AndreMenegotto/SpatialGaps

## Slope-based indicators
### Slope of species accumulation curves fitted with a Clench function: not working
````
cc1<-SACcomplet_SpName_MM(data = t500_sp) # Not working!
cc2<-SACcomplet_SEvent_MM(data = t500_sp) # Not working!
````
### Slope of species accumulation curves fitted with a methods "exact"

```
cc3<-SACcomplet_spName(data = t500_sp) # Not working!
cc4<-SACcomplet_spName(data = t500_sp) # This works
hist(cc4)
```

````
# NOT WORKING:
SACcomplet_SEvent_MM <- function(matriz)
{
  require('reshape')
  require('vegan')
  
  # Determine species unique records
  unic <- unique(matriz[,c(4,5,6,7,8)])
  
  # Create latitudinal intervals (5? bands)
  # unic[,5] <- cut(unic[,3], breaks = seq(-80, 85, 5))
  # lev <- levels(unic[,5])
  lev <-unique(matriz[,8])
  
  # Progress bar (this process may be time-consuming)
  pb <- txtProgressBar(min = 0, max = length(lev), style = 3)
  
  # Calculate level of completeness
  r <- numeric()
  
  # For each latitudinal band:
  for(i in 1:length(lev))
  {
    # Find all records
    rec <- which(unic[,5]==lev[i])
    
    if(length(rec)==0)
    {
      r[i] <- NA
    }
    else
    {
      # Create a submatrix with data of this latitudinal band
      tempMat <- unic[rec,]
      tempMat2 <- cbind(tempMat[,1:4], count=rep(1, nrow(tempMat)))
      
      # Transform the submatrix in a community matrix
      CommMat <- cast(tempMat2, latitude + longitude + date_collected ~ species_untb, value = "count")
      CommMat2 <- as.data.frame(CommMat)
      CommMat3 <- CommMat2[,-c(1,2,3)]
      CommMat3 <- as.matrix(CommMat3)
      CommMat3[which(is.na(CommMat3))] <- 0
      
      # Completeness was not calculated when there was less than 40 sampling events in the latitude because it will imply in only three points in the regression (since we use the last 10% of SACs)
      if(nrow(CommMat3) < 40)
      {
        r[i] <- NA
      }
      else
      {
        # Calculate species accumulation
        sp2 <- specaccum(CommMat3,method = "exact")
        mod1 <- try((fitspecaccum(sp2[[i]], "michaelis-menten")),
                    silent = TRUE)
        
        # estimate slope (as in Hortal et al.2004)
        a<-mod1[[10]][[1]] # Error is here
        b<-mod1[[10]][[2]] # and here
        N_rec<-length(rec)
        C<-a/(1+(b*N_rec)^2)
        
               if(is.na(C))
        {
          r[i] <- NA
        }
        else
        {
          if(C >= 0)
          {
            r[i] <- C
          }
          else
          {
            # A rapid increase of new species may create negative completeness (slope higher than 1) and we considered as zero
            r[i] <- 0
          }
        }
      }
    }
    setTxtProgressBar(pb, i)
  }  
  return(r)
}

```
### Slope of species accumulation curves fitted with a methods "exact" for "species names only": not working
