# Function to calculate mean familiarity of social neighbourhood for each focal individual

famFunction <- function(dyadData, locations){ 
  
  for(j in 1:nrow(dyadData)) {
    
    #this selects neighbours in the same grid, year, month and within 130 m of your focal individual.
    
    if (exists("neighbours")) neighbours <- NULL
    
    neighbours <- subset(locations, locations$grid==dyadData[j,grid,] & locations$year %in% dyadData[j,year,] & locations$month %in% 5 & (sqrt((30*dyadData$locX[j]-30*locations$locX)^2 + (30*dyadData$locY[j]-30*locations$locY)^2) <= 130)) # month call filters for May neighbours (neighbours at the time we're looking at)
    
    neighbours <- subset(neighbours, neighbours$squirrel_id != dyadData[j,focal,]) # remove self-self comparison
    
    neighbours <- subset(neighbours, (sqrt((30*dyadData$locX[j]-30*neighbours$locX)^2 + (30*dyadData$locY[j]-30*neighbours$locY)^2) > 0)) # remove any entries with 0 distance (as per Erin method)
    
    # now loop through neighbours to calculate individual familiarity
    
    if (exists("local")) local <- c()
    
    # including "+15" here to ensure that there is always an appropriate neighbour location to use regardless of date of specific obs. in May census
    
    fam <- subset(locations, locations$date <= (date(dyadData[j,date,])+15)) # subset of census data that occurred before May census in question
    
    if (nrow(neighbours) > 0) {
      
      for (k in 1:nrow(neighbours)) {               
        
        focalLocs <- subset(fam, fam$squirrel_id %in% dyadData[j,focal,] & (sqrt((30*dyadData$locX[j]-30*fam$locX)^2 + (30*dyadData$locY[j]-30*fam$locY)^2) <= 10)) # previous census detections of focal squirrel within 10 m of current position
        
        earliestFocal <- focalLocs[date %in% min(focalLocs$date),,] # extract earliest census observation of focal squirrel
        
        neighbourLocs <- fam[squirrel_id %in% neighbours[k,squirrel_id,]]
        
        # now need current location of neighbour
        
        neighbourLocs <- subset(neighbourLocs, (sqrt((30*neighbours$locX[k]-30*neighbourLocs$locX)^2 + (30*neighbours$locY[k]-30*neighbourLocs$locY)^2) <= 30)) # subset to locations within 10m of "current" position (relative to focal individual/year in question)
        
        # Note from Erin code: this does not take into account the possibility that the neighbouring squirrel moved from another nearby midden, but anecdotal evidence from the field suggests that anytime a neighbour moves to another nearby midden the surrounding squirrels treat that individual as a "new" neighbour)
        
        earliestNeighbour <- neighbourLocs[date %in% min(neighbourLocs$date)]
        
        firstKnownOverlap <- max(earliestFocal$date, earliestNeighbour$date)
        
        f <- as.numeric(date(dyadData[j,date,]) - firstKnownOverlap) # in days
        
        local <- c(local, f) # concatenate familiarity calculation onto temporary list that will hold all familiarity values for each focal/year combination (i.e., row in dyadData)
        
      } 
      
      set(dyadData,i = j,j = "counter", value = counter)
      
      set(dyadData,i = j,j = "numNeighbours", value = length(local))
      
      set(dyadData,i = j,j = "meanFamiliarity", value = mean(local))
      
    } else {
      
      set(dyadData,i = as.integer(j),j = "numNeighbours", value = length(local))
      
      set(dyadData,i = as.integer(j),j = "meanFamiliarity", value = NA)
      
    }
    
  }
  
  return(dyadData)
  
}


dyadData <- fullNH

dyadData <- fullNH[is.na(target),,]

j = 1

locations <- calendar_squirrelLocs



# Function to calculate familiarity of target individual in dyadic data

famFunctionTarget <- function(dyadData, locations){ 
  
  for(j in 1:nrow(dyadData)) {
    
    #this selects neighbours in the same grid, year, month and within 130 m of your focal individual.
    
    if (exists("targetID")) targetID <- NULL
    
    if (exists("target_X")) target_X <- NULL
    
    if (exists("target_Y")) target_Y <- NULL
    
    if (exists("fam")) fam <- NULL
    
    if (exists("focalLocs")) focalLocs <- NULL
    
    if (exists("earliestFocal")) earliestFocal <- NULL
    
    if (exists("targetLocs")) targetLocs <- NULL
    
    if (exists("earliestTarget")) earliestTarget <- NULL
    
    if (exists("firstKnownOverlap")) firstKnownOverlap <- NULL
    
    if (exists("f")) f <- NULL
  
    targetID <- dyadData[j,target,]
    
    target_X <- dyadData[j,t_locX,]
    
    target_Y <- dyadData[j,t_locY,]

    # including "+15" here to ensure that there is always an appropriate neighbour location to use regardless of date of specific obs. in May census
    
    fam <- subset(locations, locations$date <= (date(dyadData[j,date,])+15)) # subset of census data that occurred before May census in question (or later in that same month)
    
    # location history for focal individual
    
    focalLocs <- subset(fam, fam$squirrel_id %in% dyadData[j,focal,] & (sqrt((30*dyadData$locX[j]-30*fam$locX)^2 + (30*dyadData$locY[j]-30*fam$locY)^2) <= 10)) # previous census detections of focal squirrel within 10m of current position

    earliestFocal <- focalLocs[date %in% min(focalLocs$date),,] # extract earliest census observation of focal squirrel near same location
  
    # location history for target individual
    
    targetLocs <- fam[squirrel_id == targetID]
    
    targetLocs <- subset(targetLocs, (sqrt((30*target_X-30*targetLocs$locX)^2 + (30*target_Y-30*targetLocs$locY)^2) <= 10)) # subset to locations within 10m of "current" position (relative to focal individual/year in question)
    
    earliestTarget <- targetLocs[date %in% min(targetLocs$date)]
    
    # calculate first known overlap
    
    # Note from Erin code: this does not take into account the possibility that the neighbouring squirrel moved from another nearby midden, but evidence from the field suggests that anytime a neighbour moves to another nearby midden the surrounding squirrels treat that individual as a "new" neighbour
    
    firstKnownOverlap <- max(earliestFocal$date, earliestTarget$date) # note this requires that there is a target date!
    
    f <- as.numeric(date(dyadData[j,date,]) - firstKnownOverlap) # in days
    
    set(dyadData,i = as.integer(j),j = "targetFamiliarity", value = f)
    
  }
  
  return(dyadData)
  
}



