#Code for identifying closely related neighbours and linking their territory distance to annualized fitness measures

#Motivation: Build dataframe to investigate finding that relatedness in the social neighbourhood does not predict fitness, looking at specific relationships and continuous effects of territory distance

#Input: pairedRelations table from Script1_pedigree and various components of long-term database (dbaMidden, census, flastall, genotypes, cones)

#Start date: February 2021 by S.Walmsley

#Last modified August 2022 by S.Walmsley 


#Required packages
library(plyr)
library(dplyr)
library(magrittr)
library(lubridate)
library(evaluate)
library(data.table)
library(mefa4)
library(testthat)
library(krsp)


#For reproducibility
set.seed(9999)


# Include custom function to calculate familiarity
evaluate(file("./code/functions/famFunction.R"), keep_message = FALSE)


# Connect to database -----------------------------------

con <- krsp_connect (host = "krsp.cepb5cjvqban.us-east-2.rds.amazonaws.com",
                     
                     dbname ="krsp",
                     
                     username = Sys.getenv("krsp_user"),
                     
                     password = Sys.getenv("krsp_password"))

#secondary connection for supplementary datasets

conSupp <- krsp_connect (host = "krsp.cepb5cjvqban.us-east-2.rds.amazonaws.com",
                         
                         dbname ="krsp_suppl",
                         
                         username = Sys.getenv("krsp_user"),
                         
                         password = Sys.getenv("krsp_password"))



# Extract census records --------------------------------------------------

# Step 1: load in dbaMidden (up to 2011)

Mid <- tbl(con, "dbamidden") %>%
  
  filter(def %in% c("3","4")) %>% # including secondary middens for now
  
  collect(dbamidden)

census_1 <- data.table(Mid)

# select relevant columns for binding with new census

census_1 <- census_1[, c("reflo", "squirrel_id", "locX", "locY", "grid", "date", "def", "Sex"),]

# Step 2: load in new census (2012 onwards)

cen <- tbl(con, "census") %>%
  
  # be careful with case, database uses locx in census, but locX in dbaMidden
  
  select(reflo, squirrel_id, locx, locy, gr, census_date, sq_fate, sex) %>%
  
  # only using midden categories 1 and 2 (equivalent of extracting primary middens?)
  
  filter(sq_fate %in% c(1,2,3,4,9)) %>% # include non-midden defended areas 
  
  # not including sq_fate 6 because these secondary entries should only exist if there is a primary midden entry
  
  collect %>%
  
  mutate(locx = loc_to_numeric(locx)) # convert alphabetical locations to numeric locations

census_2 <- data.table(cen)

# rename columns to match those from dbaMidden file

colnames(census_2) <- c("reflo", "squirrel_id", "locX", "locY", "grid", "date", "sq_fate", "Sex")        

# Step 3: combine census records and clean up results

census <- rbindlist(list(census_1, census_2),use.names = TRUE,fill=TRUE)

# Clean up census file

census %<>% 
  
  mutate(date = ymd(date),
         
         year = year(date),
         
         month = month(date),
         
         locX = as.numeric(locX),
         
         locY = as.numeric(locY),
         
         reflo = as.character(reflo),
         
         Sex = as.character(Sex),
         
         grid = as.factor(grid),
         
         squirrel_id = as.character(squirrel_id)) %>%
  
  #Eliminate all entries with no reflo, LocX or LocY
  
  filter(!is.na(reflo),
         
         !is.na(locX),
         
         !is.na(locY)) %>%
  
  #Select only grids you are interested in
  
  filter(grid == "KL" | grid =="SU") %>%
  
  droplevels() %>%
  
  #Remove NA squirrel ids
  
  filter(!is.na(squirrel_id)) %>%
  
  #Exclude years outside of 2003-2014 period 
  
  filter(year > 2002 & year < 2015) %>%
  
  #Select only the columns that are relevant 
  
  select(squirrel_id, grid, date, year, month, reflo, locX, locY, Sex, def, sq_fate)


# for dbaMidden entries, identify priority territory location

census[!is.na(def),topDef:=max(def),by=c("squirrel_id","year", "month")]

census[!is.na(def),priority:=ifelse(def==topDef,"Y","N"),by=c("squirrel_id","year", "month")]

# for census entries, identify priority territory location

census[!is.na(sq_fate),topFate:=min(sq_fate),by=c("squirrel_id","year", "month")] # assuming 1 is priority if there happens to be a double entry

census[!is.na(sq_fate),priority:=ifelse(sq_fate==sq_fate,"Y","N"),by=c("squirrel_id","year", "month")]

census <- census[priority=="Y",,] 

# now that priority locations have been identified, extract one at at random for each squirrel and year

# this deals with rare instances where squirrels have two equally primary territories for the same year

squirrelLocs <- census[census[ , .I[sample(.N,1)] , by = c("squirrel_id","year", "month")]$V1]

# copy of full-calendar-year observations (best per month) to be used in familiarity calculations

calendar_squirrelLocs <- copy(squirrelLocs) 

# use May locations only for building dataframe

squirrelLocs %<>% filter(month == 5)

# remove extra dataframes that we don't require 

rm(cen,Mid,census_1,census_2)



# Calculate annual survival -----------------------------------------------

#Generate lifetime data

lifetime<-tbl(con, "flastall2") %>%
  
  # not filtering by grid, bcert, or last fate (F2) here -- looking for records of any individuals, instead using census as primary filtering step
  
  select(squirrel_id, grid=gr, sex, dates, f1, byear=byear, litter_id, dam_id, sire_id, bcert=bcert, new_grid=newgr, new_sex=newsex, datee, f2, locX, locY) %>% 
  
  collect()

#fix known errors

lifetime<-lifetime %>% 
  
  mutate (datee = ifelse(squirrel_id==6230, "1999-08-06", datee)) %>% 
  
  mutate (datee = ifelse(squirrel_id==7581, "1998-06-26", datee)) %>% 
  
  mutate (sex = ifelse (squirrel_id %in% c(6134, 5893, 5130, 4549, 4918, 7839, 5893), "M", sex)) %>% 
  
  mutate (sex = ifelse (squirrel_id==7905, "F", sex)) %>%
  
  mutate (sex = ifelse (squirrel_id==7269, "F", sex)) %>%
  
  mutate (sex = ifelse (squirrel_id==8398, "F", sex)) %>% 
  
  mutate (sex = ifelse (squirrel_id==6333, "F", sex))

lfc <- data.table(lifetime)

lfc[,dyear:=year(datee),]

lfc <- lfc[!is.na(byear),,] # remove records with unknown birth years

lfc <- lfc[!is.na(dyear),,] # remove records with unknown death years

# expand new dataset with a row for each year that each individual is alive

out <- lfc[, list(year = seq(byear, dyear), byear = byear, dyear=dyear, grid=grid, sex=sex, litter_id=litter_id, f2=f2), by = squirrel_id]

out <- out[year > 2002 & year < 2015,,] # trim down to desired study period (when we can include paternities)

out[,survival:=ifelse(year==dyear,0,1),] 

out[,age:=year-byear,] # basic age calculation

out[,unnaturalDeath:=ifelse(f2 %in% c("4", "5", "11", "12", "22") & survival==0,"Y","N"),] # IF unnatural death, remove final year of life from analysis

out <- out[unnaturalDeath=="N",,]



# Load pairwise relationship data -----------------------------------------

# Note that relationship indicates target's relation to focal ID (target --> focal)

pairedRelations <- data.table(read.csv('./input/pairwiseRelatives.csv')) # locally saved output of "pedigreeAnalyses" script

# exclude relationships that are > 0.5 relatedness but not offspring/parent/sibling (inbreeding)

pairedRelations <- pairedRelations[!is.na(relationship),,] 



# Calculate female ARS ----------------------------------------------------

evaluate(file("./code/functions/lifetimeAdapted.r"), keep_message = FALSE) # adapted portion of 'lifetime' to calculate female litter and fitness details

dfit <- data.table(dam_fit2) # convert to data table for merging 

dfit <- dfit[year > 2002 & year < 2015,,] # subset to desired years of analysis

# merge fitness values into squirrel-year dataframe

out <- merge(out, dfit, by.x = c("squirrel_id","year"), by.y = c("mother_id","year"),all.x=TRUE)

out[,squirrel_id:=as.character(squirrel_id),]

# note that female breeding as calculated here does not include "non-breeders" -- 0s are pregnancies or litters leading to no recruitment 



# Calculate male ARS ------------------------------------------------------

# bring in all siring records

flastall <- tbl(con, "flastall2") %>% 
  
  filter(gr %in% c("SU", "KL", "CH", "AG", "LL", "JO")) %>% 
  
  select(squirrel_id, gr, byear=byear, dam_id, sire_id, bcert=bcert)

flastall<-collect(flastall)

fl <- data.table(flastall)

fatheredOffspring <- fl[!is.na(sire_id),,]

hist(fatheredOffspring[,.N,by=c("byear","sire_id")]$N) 

fatheredOffspring[,maleARS:=.N,by=c("byear","sire_id")] 

fatheredOffspring = fatheredOffspring[gr %in% c("KL","SU"),,]

# sample out one record per sire per year

dads <- fatheredOffspring[fatheredOffspring[ , .I[sample(.N,1)], by=c("sire_id","byear")]$V1]

dads[,sire_id:=as.character(sire_id),]

# extract genotype table

genotypesTable <- tbl(conSupp,"genotypes") %>%
  
  collect()

genotypes <- data.table(genotypesTable)

# TEST: Do all sires actually have genotype records?

test_that("all sires have genotype records", {
  
  expect_true(all(unique(fatheredOffspring$sire_id) %in% unique(genotypes[Sex=="M",c(squirrel_id),])))
  
})

# merge fitness values into squirrel-year dataframe

out <- merge(out, dads[,c("sire_id","maleARS","byear"),], by.x = c("squirrel_id","year"), by.y = c("sire_id","byear"),all.x=TRUE)

# incorporate genotype records to attribute 0s to non-siring males

out[,genotyped:=ifelse(squirrel_id %in% genotypes$squirrel_id,"Y","N"),]

out[genotyped=="Y" & sex=="M",maleARS:=ifelse(is.na(maleARS),0,maleARS),]



# Link in territory locations ---------------------------------------------

out[,grid:=NULL,] # remove grid so we are using finer-resolution, annualized grid data from squirrelLocs

# merge annualized squirrel locations (and grid!) into squirrel-year dataframe

out <- merge(out, squirrelLocs[,c('squirrel_id','year','grid','locX','locY'),], by.x = c('squirrel_id','year'), by.y = c('squirrel_id','year'), all.x = TRUE) 



# Code mast years --------------------------------------------------------

cg <- data.table(cones_grids_years)

mastYears <- cg[Grid %in% c("KL","SU"),unique(mast),by="Year"] 

mastYears[,year:=as.integer(Year),]

mastYears[,mast:=as.factor(ifelse(V1=="n",0,1)),]

mastYears[,Year:=NULL,]

mastYears[,V1:=NULL,]

# merge mast year records into squirrel-year dataframe

out <- merge(out, mastYears, by.x="year", by.y="year", all.x=TRUE, allow.cartesian = TRUE)




# Expand pairwise relationships -------------------------------------------

final <- data.table(index = as.integer(1:20000),focal = as.character(NA), target = as.character(NA), year = as.integer(NA), relatedness = as.numeric(NA), relationship = as.character(NA))

counter = 1L

pairedRelations <- pairedRelations[focal %in% squirrelLocs$squirrel_id,,] # subset down to individuals of interest (not necessary avoids unnecessary work for the following loop)

# loop to expand possible pairwise relationships to one row per year when both individuals are alive

for (q in 1:nrow(pairedRelations)){ # for each possible pairwise close kin relationship (directional) based on pedigree
  
  for (w in out[squirrel_id==(pairedRelations[q,focal]),year,]){ # for each year that focal is alive
    
    if (w %in% out[squirrel_id==(pairedRelations[q,target,]),year,]){ # if target is also alive
      
      # add their details to the dataframe
      
      set(final,i = counter,j = "focal", value = pairedRelations[q,focal,])
      
      set(final,i = counter,j = "target", value = pairedRelations[q,target,])
      
      set(final,i = counter,j = "year", value = w)
      
      set(final,i = counter,j = "relatedness", value = pairedRelations[q,relatedness,])
      
      set(final,i = counter,j = "relationship", value = pairedRelations[q,relationship,])
      
      counter = as.integer(counter + 1)
      
    }
    
  }
  
}

final <- final[!is.na(focal),,] # remove empty rows in dataframe

final[,dyad:=paste(sort(c(focal,target))[1],sort(c(focal,target))[2],sep="_"),by=index] # alphabetical dyad name (works bi-directionally)



# Merge in focal data -----------------------------------------------------

final_2 <- merge(final, out, by.x=c("focal","year"), by.y=c("squirrel_id","year"),all.x=TRUE)

final_2 <- final_2[!is.na(locX),,] # exclude entries where squirrel must be alive but does not have a recorded location  

# Copy dataset and rename columns so "target" individual details can be incorporated smoothly

targetData <- copy(out[,c("squirrel_id","locX","locY","sex","year","grid","litter_id","age")])

colnames(targetData) <- c("squirrel_id", "t_locX", "t_locY", "t_sex", "year","t_grid","t_litter_id","t_age")

# merge target details into main dataset

final_2 <- merge(final_2,targetData, by.x=c("target","year"),by.y=c("squirrel_id","year"))



# Calculate distance between territories ----------------------------------

final_2[,distance:= sqrt((30*(locX)-30*(t_locX))^2 + (30*(locY)-30*(t_locY))^2),]

# exclude any records when the dyad spans multiple grids

distances <- final_2[!is.na(distance) & grid==t_grid,,] 

# littermates

distances[,is_mother:=ifelse(relationship=="parent" & t_sex=="F","Y","N"),]

distances[is_mother=="Y",littermates:="Y",by="index"]

distances[is_mother=="N",littermates:=ifelse(litter_id==t_litter_id,"Y","N"),by="index"]

# fine-grain relationships

distances[,relationshipFine:=ifelse(relationship=="parent" & t_sex=="F","mother",NA),]

distances[is.na(relationshipFine),relationshipFine:=ifelse(relationship=="parent" & t_sex=="M","father",NA),]

distances[is.na(relationshipFine),relationshipFine:=ifelse(relationship=="offspring" & t_sex=="F","daughter",NA),]

distances[is.na(relationshipFine),relationshipFine:=ifelse(relationship=="offspring" & t_sex=="M","son",NA),]

distances[is.na(relationshipFine) & relationship=="sibling",relationshipFine:="sibling",]

distances[sex=="F",.N,by=is.na(all_litters_fit)] 

distances[sex=="M",.N,by=is.na(maleARS)] 

## Subset "nearest" dataset based on each focal individual's nearest close kin, by year

distances[,minDistance:=min(distance),by=c("focal","year")]

nearest <- distances[distance==minDistance,,] 

nearest[,tieBreak:=1:.N,by=c("focal",'year')]

# select just one of these records

nearest <- nearest[tieBreak==1,,]

# add a column of log-transformed linear distance

nearest[,logDist:=log10(distance),]

# convert relationships to factor

nearest[,relationship:=as.factor(relationship),]

nearest[,relationshipFine:=as.factor(relationshipFine),]


# TEST: All records in "nearest" dataframe have distances and known relationship?

test_that("Nearest records are complete", {
  
  expect_true(nearest[is.na(relationship),.N,]==0)
  
  expect_true(nearest[is.na(relationshipFine),.N,]==0)
  
  expect_true(nearest[is.na(locX),.N,]==0)
  
  expect_true(nearest[is.na(t_locX),.N,]==0)
  
  expect_true(nearest[is.na(distance),.N,]==0)
  
})



# Incorporate "no-kin" squirrels ------------------------------------------

outLocs <- out[!is.na(locX),,] # basically squirrelLocs with fitness measures added

lifeStatus <- merge(outLocs, nearest, by.x=c("squirrel_id","year"),by.y=c("focal","year"),all.x=TRUE)

# 732 squirrels (just 425 for nearest)

lifeStatus[,closeKinPresent:=ifelse(is.na(t_locX),"N","Y"),] # adding 816 observations

lifeStatus[,.N,by="closeKinPresent"]

# 838 observations where no close kin are present





# No-kin dataset: social neighbourhood ------------------------------------


#these individuals will be listed as "unrelated" (meaning not close kin), with distance to a random individual within 130m 

noKin <- lifeStatus[closeKinPresent=="N",,]

# Calculating neighbourhoods

noKin[,grid:=grid.x,]

noKin[,locX:=locX.x,]

noKin[,locY:=locY.x,]

noKin[,neighbourDistance:=as.numeric(NA),]

noKin[,neighbourID:=as.character(NA),]

#benefiting from adaptation of Erin's social neighbourhood code

for (j in 1:nrow(noKin)) { # looping through rows of 'non-kin' individuals
  
  #this selects neighbours in the same grid, year, month and within 130 m of your focal individual.
  
  if (exists("neighbours")) neighbours <- NULL
  
  neighbours <- subset(noKin, noKin$grid==noKin$grid[j] & noKin$year==noKin$year[j] & (sqrt((30*noKin$locX[j]-30*noKin$locX)^2 + (30*noKin$locY[j]-30*noKin$locY)^2) <= 130))
  
  neighbours <- subset(neighbours, neighbours$squirrel_id != noKin[j,squirrel_id,]) # remove self-self comparison!
  
  if (nrow(neighbours)>=1) { # if there is a neighbour
    
    index <- sample(1:nrow(neighbours),1) # sample a neighbour at random
    
    set(noKin,i = j,j = "neighbourDistance", value = (sqrt((30*noKin[j,locX,]-30*neighbours[index,locX,])^2 + (30*noKin[j,locY,]-30*neighbours[index,locY,])^2)))
    
    set(noKin,i = j,j = "neighbourID", value = neighbours[index,squirrel_id,])
    
    set(noKin,i = j, j = "t_locX", value = neighbours[index,locX,]) # added July 2022 for familiarity calculations
    
    set(noKin,i = j, j = "t_locY", value = neighbours[index,locY,]) # added July 2022 for familiarity calculations
    
  } else { 
    
    set(noKin,i = j,j = "neighbourDistance", value = NA)
    
    set(noKin,i = j,j = "neighbourID", value = NA)
    
    set(noKin,i = j, j = "t_locX", value = NA) # added July 2022 for familiarity calculations
    
    set(noKin,i = j, j = "t_locY", value = NA) # added July 2022 for familiarity calculations
    
  }
  
}

# check out distribution

hist(noKin$neighbourDistance)

#nK <- noKin[,c('squirrel_id','neighbourID','grid.x','neighbourDistance','year','sex.x','age.x','mast.x','survival.x','all_litters_fit.x','maleARS.x','relationship')]

nK <- noKin[,c('squirrel_id','neighbourID','grid.x','neighbourDistance','year','sex.x','age.x','mast.x','survival.x','all_litters_fit.x','maleARS.x','relationship','locX','locY','t_locX','t_locY')]

colnames(nK) <- c('focal','target','grid','distance','year','sex','age','mast','survival','all_litters_fit','maleARS','relationshipFine', 'locX', 'locY', 't_locX', 't_locY')

nK[,relationshipFine:="unrelated",]

nK[,logDist:=log10(distance),]

nK[,ii:=.I,]

nK[,dyad:=paste(sort(c(focal,target))[1],sort(c(focal,target))[2],sep="_"),by=ii] # alphabetical dyad name (works bi-directionally)

nHood <- nearest[distance<=130,,]

fullNH <- rbindlist(list(nHood,nK),fill=TRUE)

fullNH$relationshipFine <- factor(fullNH$relationshipFine,levels=c('unrelated','daughter','mother', 'sibling', 'son','father'))

levels(fullNH$relationshipFine) <- c("Non-kin", "Daughter", "Mother", "Sibling", "Son", "Father")


# Add mean familiarity to neighbourhood dataset ---------------------------

fullNH[,date:=paste(year,"-05-15",sep="")] # add date column for May census

fullNH[,numNeighbours:=as.integer(NA),]

fullNH[,meanFamiliarity:=as.numeric(NA),]

#fullNH <- famFunction(fullNH, calendar_squirrelLocs)

# Add target familiarity to neighbourhood datasets ------------------------

fullNH[,targetFamiliarity:=as.numeric(NA),]

fullNH <- famFunctionTarget(fullNH, calendar_squirrelLocs)

fullNH[targetFamiliarity<0,targetFamiliarity:=0,] # fix any negative familiarity values (same as Erin code, happens if census dates are different days)

# Save output as .csv -----------------------------------------------------

#write.csv(fullNH,"./output/socialNeighbourhood_August2022.csv")





# No-kin dataset: all distances -------------------------------------------


#these individuals will be listed as "unrelated"(meaning not close kin), with distance to a random individual within 130m 

noKin <- lifeStatus[closeKinPresent=="N",,]

# Calculating neighbourhoods

noKin[,grid:=grid.x,]

noKin[,locX:=locX.x,]

noKin[,locY:=locY.x,]

noKin[,neighbourDistance:=as.numeric(NA),]

noKin[,neighbourID:=as.character(NA),]

#benefiting from adaptation of Erin's social neighbourhood code

for (j in 1:nrow(noKin)) {
  
  #this selects neighbours in the same grid, year, month and within 130 m of your focal individual.
  
  if (exists("neighbours")) neighbours <- NULL
  
  neighbours <- subset(noKin, noKin$grid==noKin$grid[j] & noKin$year==noKin$year[j] & (sqrt((30*noKin$locX[j]-30*noKin$locX)^2 + (30*noKin$locY[j]-30*noKin$locY)^2) <= 10000))
  
  neighbours <- subset(neighbours, neighbours$squirrel_id != noKin[j,squirrel_id,]) # remove self-self comparison!
  
  if (nrow(neighbours)>=1) {
    
    index <- sample(1:nrow(neighbours),1) # sample a neighbour at random
    
    #noKin[j,neighbourDistance:=(sqrt((30*noKin$locX[j]-30*neighbours$locX[index])^2 + (30*noKin$locY[j]-30*neighbours$locY[index])^2)),]
    
    set(noKin,i = j,j = "neighbourDistance", value = (sqrt((30*noKin[j,locX,]-30*neighbours[index,locX,])^2 + (30*noKin[j,locY,]-30*neighbours[index,locY,])^2)))
    
    set(noKin,i = j,j = "neighbourID", value = neighbours[index,squirrel_id,])
    
    set(noKin,i = j, j = "t_locX", value = neighbours[index,locX,]) # added July 2022 for familiarity calculations
    
    set(noKin,i = j, j = "t_locY", value = neighbours[index,locY,]) # added July 2022 for familiarity calculations
    
  } else {
    
    set(noKin,i = j,j = "neighbourDistance", value = NA)
    
    set(noKin,i = j,j = "neighbourID", value = NA)
    
    set(noKin,i = j, j = "t_locX", value = NA) # added July 2022 for familiarity calculations
    
    set(noKin,i = j, j = "t_locY", value = NA) # added July 2022 for familiarity calculations
    
  }
  
}

# check out distribution

hist(noKin$neighbourDistance)

noKin[,.N,by=is.na(neighbourDistance)] 

noKin[neighbourDistance==0,.N,]

nK <- noKin[,c('squirrel_id','neighbourID','grid.x','neighbourDistance','year','sex.x','age.x','mast.x','survival.x','all_litters_fit.x','maleARS.x','relationship','locX','locY','t_locX','t_locY')]

colnames(nK) <- c('focal','target','grid','distance','year','sex','age','mast','survival','all_litters_fit','maleARS','relationshipFine','locX','locY','t_locX','t_locY')

nK[,relationshipFine:="unrelated",]

nK[,logDist:=log10(distance),]

nK[,ii:=.I,]

nK[,dyad:=paste(sort(c(focal,target))[1],sort(c(focal,target))[2],sep="_"),by=ii] # alphabetical dyad name (works bi-directionally)

fullAll <- rbindlist(list(nearest,nK),fill=TRUE)

fullAll$relationshipFine <- factor(fullAll$relationshipFine,levels=c('unrelated','daughter','mother', 'sibling', 'son','father'))

levels(fullAll$relationshipFine) <- c("Non-kin", "Daughter", "Mother", "Sibling", "Son", "Father")


# Add mean familiarity to entire study area dataset -----------------------

fullAll[,date:=paste(year,"-05-15",sep="")] # add date column for May census

fullAll[,numNeighbours:=as.integer(NA),]

fullAll[,meanFamiliarity:=as.numeric(NA),]

fullAll <- famFunction(fullAll, calendar_squirrelLocs)


# Add target familiarity to entire study area datasets ------------------------

fullAll[,targetFamiliarity:=as.numeric(NA),]

fullAll <- famFunctionTarget(fullAll, calendar_squirrelLocs)

fullAll[targetFamiliarity<0,targetFamiliarity:=0,] # fix any negative familiarity values (same as Erin code, happens if census dates are different days)

# Save output as .csv -----------------------------------------------------

#write.csv(fullAll,"./output/entireGrid_August2022.csv")


# Note that all columns unnecessary for subsequent analysis (Script 3) were removed prior to creation of final dataset


