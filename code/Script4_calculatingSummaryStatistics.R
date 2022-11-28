#Code for estimating influence of closely related conspecifics on fitness in red squirrels

#Motivation: Investigate finding that relatedness in the social neighbourhood does not predict fitness, looking at specific relationships and continuous effects of territory distance

#Input: dataframes built and tested in Script2_dataframe 

#Start date: March 2021 by S.Walmsley

#Last modified July 2022 by S.Walmsley 


#Required packages
library(plyr)
library(dplyr)
library(magrittr)
library(lubridate)
library(evaluate)
library(data.table)
library(mefa4)
library(ggplot2)
library(testthat)
library(paletteer)
library(glmmTMB)



# Import saved data -------------------------------------------------------


allDistances <- data.table(read.csv('./input/kinDistances_All_ForSummaryStats.csv'))

allDistances$relationshipFine <- factor(allDistances$relationshipFine,levels=c("Non-kin", "daughter", "mother", "sibling", "son", "father"))

levels(allDistances$relationshipFine) <- c("'Non-kin'", "Daughter", "Mother", "Sibling", "Son", "Father")


#fullNH <- data.table(read.csv('./output/socialNeighbourhood_July2022.csv'))

fullNH <- data.table(read.csv('./input/socialNeighbourhood_final.csv'))

fullNH$relationshipFine <- factor(fullNH$relationshipFine,levels=c("Non-kin", "Daughter", "Mother", "Sibling", "Son", "Father"))

levels(fullNH$relationshipFine) <- c("'Non-kin'", "Daughter", "Mother", "Sibling", "Son", "Father")


#fullGrid <- data.table(read.csv('./output/entireGrid_August2022.csv'))

#fullGrid <- data.table(read.csv('./output/entireGrid_July2022_10k.csv'))

#fullGrid <- data.table(read.csv('./output/entireGrid_July2022.csv'))

fullGrid <- data.table(read.csv('./input/entireGrid_final.csv'))

fullGrid$relationshipFine <- factor(fullGrid$relationshipFine,levels=c("Non-kin", "Daughter", "Mother", "Sibling", "Son", "Father"))

levels(fullGrid$relationshipFine) <- c("'Non-kin'", "Daughter", "Mother", "Sibling", "Son", "Father")

fgNK <- fullGrid[relationshipFine=="'Non-kin'",,]




### CODE OPTION A ###

## Operations social neighbourhood =================================

allDistances <- allDistances[distance <= 130,,]

#allDistances[,numYears:=length(unique(year)),by="dyad"] # calculate number of years we have records for each dyad

#uniqueDyads <- allDistances[allDistances[ , .I[sample(.N,1)], by=c("dyad")]$V1]

#uniqueDyads[,mean(numYears),] # 1.407

#uniqueDyads[,range(numYears)] # 1 4

adNH <- allDistances[distance<=130,,]

hist(adNH$distance)

mean(adNH[,length(unique(year)),by="dyad"]$V1) # 1.407

range(adNH[,length(unique(year)),by="dyad"]$V1) # 1 4


aD <- rbindlist(list(adNH,fgNK),fill = TRUE)

hist(aD$distance)

#mean(aD[,length(unique(year)),by="dyad"]$V1) # mean 1.113

#range(aD[,length(unique(year)),by="dyad"]$V1) # 1 4


aD[relationshipFine=="'Non-kin'",numKin:=0,by=c('focal','year')]

aD[relationshipFine!="'Non-kin'",numKin:=length(unique(target)),by=c('focal','year')]

mean(aD[,numKin,by=c('focal','year')]$numKin) # 0.9302615

median(aD[,numKin,by=c('focal','year')]$numKin) # 0

range(aD$numKin) # 0 6

aD[,dd:=ifelse(numKin>0,1,0),by=c('focal','year')] # 1 if have kin


# need to narrow in on one record per individual per year

uA <- aD[aD[ , .I[sample(.N,1)] , by = c("focal","year")]$V1]

mean(uA$numKin) ############################################### mean 0.564291

range(uA$numKin) ############################################## range 0-6

hist(uA$numKin) 



### CODE OPTION B ###

## Operations full grid =================================


allDistances[,numYears:=length(unique(year)),by="dyad"] # calculate number of years we have records for each dyad

# fix this to use squirrelLocs, with one record per dyad

uniqueDyads <- allDistances[allDistances[ , .I[sample(.N,1)], by=c("dyad")]$V1]

uniqueDyads[,mean(numYears),] # 1.474273

uniqueDyads[,range(numYears)] # 1 4


aD <- rbindlist(list(allDistances,fgNK),fill = TRUE)

hist(aD$distance)

mean(aD[,length(unique(year)),by="dyad"]$V1) # mean 1.17806

range(aD[,length(unique(year)),by="dyad"]$V1) # range 1-4


aD[relationshipFine=="'Non-kin'",numKin:=0,by=c('focal','year')]

aD[relationshipFine!="'Non-kin'",numKin:=length(unique(target)),by=c('focal','year')]

mean(aD[,numKin,by=c('focal','year')]$numKin) # 1.44898

median(aD[,numKin,by=c('focal','year')]$numKin) # 0  # 1

range(aD$numKin) # 0 7

aD[,dd:=ifelse(numKin>0,1,0),by=c('focal','year')] # 1 if have kin


# need to narrow in on one record per individual per year

uA <- aD[aD[ , .I[sample(.N,1)] , by = c("focal","year")]$V1]

mean(uA$numKin) ############################################### mean 0.8294525

range(uA$numKin) ############################################## range 0-6 # 0 7

hist(uA$numKin) 




# sex distance model

#dt <- allDistances[allDistances[,.I[sample(.N,1)],by=c("dyad")]$V1]

dt <- fullGrid[relationshipFine != "'Non-kin'",,]

dt <- dt[dt[,.I[sample(.N,1)],by=c("dyad")]$V1]

# note that the specific coefficients of the following model will vary slightly depending on which observations are sampled

m <- glmmTMB(data=dt[distance>0],distance~sex+(1|focal),family=Gamma(link="log"))

summary(m)

mdf <- ggpredict(m, terms = c("sex"))

mdf







###########################################################################

# Boxplot script ----------------------------------------------------------


std <- function(x) sd(x)/sqrt(length(x))


allDistances[,sexRF:=as.factor(paste(relationshipFine,sex)),]

allDistances$sexRF <- factor(allDistances$sexRF,levels=c("Daughter F", "Mother F", "Sibling F", "Son F", "Father F","Daughter M", "Mother M", "Sibling M", "Son M", "Father M"))



allDistances[,sexMed:=median(distance),by="sex"]

allDistances[,rfMed:=median(distance),by="sexRF"]

allDistances[,printSex:=ifelse(sex=="F","Female perspective","Male perspective"),]



allDistances[,median(distance),by="sex"]

allDistances[,sd(distance),by="sex"]

allDistances[,std(distance),by="sex"]


sub <- allDistances[allDistances[ , .I[sample(.N,1)], by=c("focal","target")]$V1]

sub[,sexMed:=median(distance),by="sex"]

sub[sex=="F",median(distance),by=c("relationshipFine")]


golds <- c('#FFDC73','#FFCC00','#ECBD00','#BF9B30','#B8860B')

fcols <- c('grey60',paletteer_dynamic('cartography::green.pal',6)[2:6])

mcols <- c('grey60',golds)

#acols <- c(fcols[2:6],mcols[2:6])

acols <- c(fcols[2],mcols[4])



gx <- ggplot(data=sub,aes(x=as.factor(relationshipFine),y=distance,fill=printSex,color=printSex))+
  
  facet_wrap(~printSex)+
  
  geom_hline(aes(yintercept = sexMed), color = "grey50", size =1.5) +
  
  geom_point(aes(color=printSex),alpha=0.5,position = position_jitter(width = 0.15,height = 0), size = 2, shape=21)+
  
  geom_boxplot(outlier.shape = NA, alpha = 0.0, width = 0.4, colour = "black",size=0.9,notch=TRUE)+
  
  #scale_color_manual(shape = NA, alpha = 0.0, width = 0.4, size=0.725,colour = "black",notch=TRUE)+
  
  #geom_point(aes(y = rfMed),color="black", size = 5)+
  
  #geom_flat_violin(position = position_nudge(x = .25, y = 0), adjust = 1, trim = TRUE, alpha = 1, colour = NA)+
  
  guides(fill = "none",colour="none")+
  
  xlab('')+ylab('Distance to territories of primary kin (m)')+
  
  scale_fill_manual(values=acols)+
  
  #ylim(0,200)+
  
  #geom_segment(aes(xend=relationshipFine,y=sexMed,yend=rfMed),color="black",size=1.5)+
  
  #geom_segment(aes(x = sexRF, xend = sexRF,y = sexMed, yend = rfMed),size = 0.8) +
  
  theme(plot.title = element_text(size = 16),strip.background = element_rect(fill="white"),axis.ticks = element_blank(),axis.text=element_text(size=12,hjust=0.35),panel.background=element_blank(),panel.border=element_rect(fill=NA,colour = 'grey50'),panel.grid.major.y = element_line(colour = "grey95",size = 1),panel.grid.minor.y = element_line(colour = "grey95",size = 1),strip.text=element_text(size="12"))

gx

#png("./figures/FigX_boxplot_Nov2022.png",width=8,height=6,units='in',res=700)
gx
#dev.off()





