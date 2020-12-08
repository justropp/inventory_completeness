# Indicators of inventory completeness for global biodiversity databases
Juliana Stropp<br/>

This document presents the R code used to calculate indicators of inventory completeness presented in Stropp et al. (in prep.; Table 1).<br/>

R functions were adapted from: https://github.com/AndreMenegotto/SpatialGaps

## Slope-based indicators
#### 1. Slope of species accumulation curves fitted with a Clench function
````
cc1<-SACcomplet_SpName_MM(matriz = t500_sp) # Not working!
cc2<-SACcomplet_SEvent_MM(matriz = t500_sp) # Not working!
````
#### 2. Slope of species accumulation curves fitted the method "exact"

```
cc3<-SACcomplet_spName(matriz = t500_sp) # Not working!
cc4<-SACcomplet_SEvent(matriz = t500_sp) # This works
hist(cc4)
```
## Abundance-based indicators
#### 1. Chao2

```
cc5<-Chao2_spName(matriz = t500_sp) # Not working!
cc6<-Chao2_SampEvent(matriz = t500_sp) # This works
```
#### 2. Chao2, with its adjustment for small sample sizes
```
cc7<-Chao2adj_spName(matriz = t500_sp) # Not working!
cc8<-Chao2adj_SEvent(matriz = t500_spp) # This works
```
#### 3. Chao1
```
cc9<-Chao1_spName (t500_spp) # Not working
cc10<-Chao1_SEvent (t500_spp) # This works
```
#### 4. Sample coverage (Chao and Jost, 2012)
````
cc11<-ChaoJost_spName(t500_spp) #Not working
cc12<-ChaoJost_SEvent(t500_spp) # This works
````
```
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
Slope - method "exact"; using species name only
````
# NOT WORKING

````

Slope - method "exact"; using sampling events

````
# This works
SACcomplet_SEvent <- function(matriz)
{
  require('reshape')
  require('vegan')
  
  # Select unique observation for each species for each group (cell_id and sampling_intensity)
  unic <- unique(matriz[,c(4,5,6,7,8)])
  
  # Define groups to calculate completeness: cell_id, synth_com, sampling_intensity
  lev <-unique(matriz[,8])
  
  # Progress bar (this process may be time-consuming)
  pb <- txtProgressBar(min = 0, max = length(lev), style = 3)
  
  # Start to calculate inventory completeness
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
      # Create a submatrix with data in each group
      tempMat <- unic[rec,]
      tempMat2 <- cbind(tempMat[,1:4], count=rep(1, nrow(tempMat)))
      
      # Transform the submatrix in a community matrix
      CommMat <- cast(tempMat2, latitude + longitude + date_collected ~ species_untb, value = "count")
      CommMat2 <- as.data.frame(CommMat)
      CommMat3 <- CommMat2[,-c(1,2,3)]
      CommMat3 <- as.matrix(CommMat3)
      CommMat3[which(is.na(CommMat3))] <- 0
      
      # Completeness was not calculated when less than 40 sampling events 
      # is recorded in a cell. For such cases only three points would be used 
      # in the regression to extract mean slope of the last 10% of SACs
      if(nrow(CommMat3) < 40)
      {
        r[i] <- NA
      }
      else
      {
        # Calculate species accumulation
        sp1 <- specaccum(CommMat3, method = "exact")
        
        # Extract the last 10%
        x <- quantile(sp1$richness, probs = c(0.9))
        Lten <- which(sp1$richness >= x)
        vec <- sp1$richness[Lten]
        
        # Result (1 - regression slope)
        slope <- 1-lm(vec~seq(1,length(vec),1))[[1]][[2]]
        
        if(is.na(slope))
        {
          r[i] <- NA
        }
        else
        {
          if(slope >= 0)
          {
            r[i] <- slope
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
````

