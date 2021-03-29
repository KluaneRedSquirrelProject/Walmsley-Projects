#Code for expanding pedigree into pairwise relationships of close relatives

#Motivation: Will produce output file for use in 'closeRelatives'

#Start date: February 11 2021 by S.Walmsley

#Last modified March 15 2021 by S.Walmsley 

#Required packages
library(plyr)
library(dplyr)
library(data.table)
library(krsp)
library(kinship2)
library(devtools)
library(pedantics)
library(Matrix)

#For reproducibility
set.seed(1234)

## Establish connection to KRSP database =================================

con <- krsp_connect (host = "krsp.cepb5cjvqban.us-east-2.rds.amazonaws.com",
                     dbname ="krsp",
                     username = Sys.getenv("krsp_user"),
                     password = Sys.getenv("krsp_password"))
#run create_pedigree script from github
source_url("https://raw.githubusercontent.com/KluaneRedSquirrelProject/krsp-functions/master/create_pedigree.R")


## Calculate pairwise relatedness matrix =================================

# run pedantics function on pedigree
krspPedigreeSummary<-pedigreeStats(krsp_pedigree, graphicalReport="n")

# extract relatedness matrix
aa <- krspPedigreeSummary$Amatrix# 400 million elements, mostly 0s
# convert into sparse matrix for easier processsing
aaSparse <- Matrix(aa,sparse=TRUE)


## Expand relatedness matrix into long dyadic dataframe =================================

# convert pedigree to data table for fast searching 
dt_krsp_pedigree <- data.table(krsp_pedigree)

#use function for melting sparse matrices (verbose matrix for entire pedigree is too large to melt -- 400 million + elements)
dt <- data.table(setNames(mefa4::Melt(aaSparse), c('focal', 'target', 'relatedness'))) # this works

#subset to create new data table including only close relatives and excluding self-self comparisons
relatives <- dt[relatedness>=0.5 & focal!=target,,]
relatives[,index:=.I,] # add generic index column, useful for next steps


## Assign basic genealogical relationships =================================

#identify and label parent relationships [this takes a little while ~20 minutes on my machine]
relatives[,relationship:=ifelse(target %in% dt_krsp_pedigree[id==focal,c(as.character(dam),as.character(sire)),],"parent",NA),by=index] 

#identify and label offspring relationships
relatives[is.na(relationship),relationship:=ifelse(focal %in% dt_krsp_pedigree[id==target,c(as.character(dam),as.character(sire)),],"offspring",NA),by=index] 

#identify and label full sibling relationships
relatives[is.na(relationship),relationship:=ifelse(all(dt_krsp_pedigree[id==focal,c(as.character(dam),as.character(sire)),] %in% dt_krsp_pedigree[id==target,c(as.character(dam),as.character(sire)),]),"sibling",NA),by=index] 

# can now write out .csv file for saving or simply use dataframe in other scripts

