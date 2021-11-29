#Code for estimating influence of closely related conspecifics on fitness in red squirrels

#Motivation: Investigate finding that relatedness in the social neighbourhood does not predict fitness, looking at specific relationships and continuous effects of territory distance

#Input: dataframes built and tested in Script2_dataframe 

#Start date: March 10 2021 by S.Walmsley

#Last modified November 22 2021 by S.Walmsley 

#Required packages
library(glmmTMB)
library(ggeffects)
library(mgcv)
library(plyr)
library(dplyr)
library(magrittr)
library(lubridate)
library(evaluate)
library(data.table)
library(mefa4)
library(testthat)
library(krsp)
library(ggplot2)
library(paletteer)
library(DHARMa)
library(ggsci)
library(viridis)
library(ggtext)


#Custom function to print outputs from models for easy pasting into word
evaluate(file("./functions/printOutput.R"), keep_message = FALSE) # adapted portion of 'lifetime' to calculate female litter and fitness details
#Custom function to create plots using ggplot2
evaluate(file("./functions/plotGAMs.R"), keep_message = FALSE) # adapted portion of 'lifetime' to calculate female litter and fitness details
#GeomFlatViolin available at https://rdrr.io/github/jorvlan/openvis/src/R/R_rainclouds.R
evaluate(file("./functions/GeomFlatViolin.R"), keep_message = FALSE) # adapted portion of 'lifetime' to calculate female litter and fitness details


#Colors for figures
fcols <- c('grey60',"#247524FF","#247524FF","#247524FF","#247524FF","#247524FF")
mcols <- c("grey60","#B8860B","#B8860B","#B8860B","#B8860B" ,"#B8860B" )


# Import saved data -------------------------------------------------------

fullNH <- data.table(read.csv('./output/nonKinDataset_socialNeighbourhood2.csv'))
fullNH$relationshipFine <- factor(fullNH$relationshipFine,levels=c("Non-kin", "Daughter", "Mother", "Sibling", "Son", "Father"))
levels(fullNH$relationshipFine) <- c("'Non-kin'", "Daughter", "Mother", "Sibling", "Son", "Father")

fullGrid <- data.table(read.csv('./output/nonKinDataset_allDistances2.csv'))
fullGrid$relationshipFine <- factor(fullGrid$relationshipFine,levels=c("Non-kin", "Daughter", "Mother", "Sibling", "Son", "Father"))
levels(fullGrid$relationshipFine) <- c("'Non-kin'", "Daughter", "Mother", "Sibling", "Son", "Father")


# Female survival ---------------------------------------------------------

# color for visualization
vcols <- fcols
# calculate sample sizes
fgNo0 <- fullGrid[distance>0,,]
fgNo0[,numByRelation:=.N,by=c("sex","relationshipFine")]
fgNo0[,group_col:=relationshipFine,]
sampleSizes <- fgNo0[fgNo0[ , .I[sample(.N,1)] , by = c("sex","group_col")]$V1]

# non-linear visualization for full grid dataset
m <- gam(data=fullGrid[sex=="F" & distance>0],survival ~ relationshipFine + s(logDist,by=relationshipFine,k=3)+age+I(age^2)+mast,binomial)
summary.gam(m)
#png("./figures/Nov_femaleSurvivalGAMoutput.png",width=9,height=7,units='in',res=700)
plotGAM_viridis(m,sampleSizes,"F","Influence of nearby kin on female survival","Probability of survival",1,2.6,0.08)
#dev.off()

# social neighborhood model
m1 <- glmmTMB(data=fullNH[sex=="F" & distance>0],survival~logDist*relationshipFine+age+I(age^2)+as.factor(mast)+(1|focal)+(1|dyad)+(1|year),family="binomial")
summary(m1)
plot(ggpredict(m1,terms = c("logDist","relationshipFine")),facet=TRUE,add.data=TRUE)
plot(simulateResiduals(m1, plot = F)) # basic validation check
#printOutput(m1,"Nov_femaleSurvival")


# full grid model
m2 <- glmmTMB(data=fullGrid[sex=="F" & distance>0],survival~relationshipFine*poly(logDist,2)+age+I(age^2)+mast+(1|focal)+(1|dyad)+(1|year),family="binomial")
summary(m2)
plot(simulateResiduals(m2, plot = F)) # basic validation check
#printOutput(m2,"Nov_femaleSurvival_fullGrid")

# calculate magnitude of daughter effect
mdf <- ggpredict(m1, terms = c("logDist [1,1.477121,2,2.113943]", "relationshipFine [Daughter]", "mast [0]", "age [2]"))
mdf

# Female reproductive success ---------------------------------------------

# color for visualization
vcols <- fcols
# calculate sample sizes
fgNo0 <- fullGrid[distance>0 & !is.na(all_litters_fit),,]
fgNo0[,numByRelation:=.N,by=c("sex","relationshipFine")]
fgNo0[,group_col:=relationshipFine,]
sampleSizes <- fgNo0[fgNo0[ , .I[sample(.N,1)] , by = c("sex","group_col")]$V1]

# non-linear visualization for full grid dataset
vcols <- fcols
m <- gam(data=fullGrid[sex=="F" & distance>0],all_litters_fit ~ relationshipFine + s(logDist,by=relationshipFine,k=3)+age+I(age^2)+mast,poisson)
summary.gam(m)
#png("./figures/Nov_femaleARSGAMoutput.png",width=9,height=7,units='in',res=700)
plotGAM_viridis(m,sampleSizes,"F","Influence of nearby kin on female reproductive success","Pups recruited",7,1.38,6.8)
#dev.off()

# Sensitivity analysis for inclusion of non-breeding females

# calculate sample sizes
fgNo0 <- fullGrid[distance>0,,]
fgNo0[,numByRelation:=.N,by=c("sex","relationshipFine")]
fgNo0[,group_col:=relationshipFine,]
sampleSizes <- fgNo0[fgNo0[ , .I[sample(.N,1)] , by = c("sex","group_col")]$V1]

fullGrid[sex=="F",nbARS:=ifelse(is.na(all_litters_fit),0,all_litters_fit),]
# non-linear visualization for full grid dataset including non-breeders as 0s
mS <- gam(data=fullGrid[sex=="F" & distance>0],nbARS ~ relationshipFine + s(logDist,by=relationshipFine,k=3)+age+I(age^2)+mast,poisson)
summary.gam(mS)
#png("./figures/Nov_breeders_femaleARSGAMoutput.png",width=9,height=7,units='in',res=700)
plotGAM_viridis(mS,sampleSizes,"F","Influence of nearby kin on female reproductive success","Pups recruited",7,1.38,6.8)
#dev.off()

# social neighborhood model
m1 <- glmmTMB(data=fullNH[sex=="F" & distance>0],all_litters_fit~logDist*relationshipFine+age+I(age^2)+mast+(1|focal) +(1|dyad) + (1|year),family="poisson")
summary(m1)
plot(simulateResiduals(m1, plot = F)) # basic validation check
#printOutput(m1,"Nov_femaleARS")

# full grid model
# no apparent non-linearities in full grid dataset -- model linearly
m2 <- glmmTMB(data=fullGrid[sex=="F" & distance>0],all_litters_fit~logDist*relationshipFine+age+I(age^2)+mast+(1|focal) +(1|dyad)+ (1|year),family="poisson")
summary(m2)
plot(simulateResiduals(m2, plot = F)) # basic validation check
#printOutput(m2,"Nov_femaleARS_fullGrid")



# Male survival -----------------------------------------------------------

# color for visualization
vcols <- mcols
# calculate sample sizes
fgNo0 <- fullGrid[distance>0,,]
fgNo0[,numByRelation:=.N,by=c("sex","relationshipFine")]
fgNo0[,group_col:=relationshipFine,]
sampleSizes <- fgNo0[fgNo0[ , .I[sample(.N,1)] , by = c("sex","group_col")]$V1]

# non-linear visualization for full grid dataset
m <- gam(data=fullGrid[sex=="M" & distance>0],survival ~ relationshipFine + s(logDist,by=relationshipFine,k=3)+age+I(age^2)+mast,binomial)
summary.gam(m)
#png("./figures/Nov_maleSurvivalGAMoutput.png",width=9,height=7,units='in',res=700)
plotGAM_viridis(m,sampleSizes,"M","Influence of nearby kin on male survival","Probability of survival",1,1.55,0.08)
#dev.off()

# social neighborhood model
m1 <- glmmTMB(data=fullNH[sex=="M" & distance>0],survival~logDist*relationshipFine+age+I(age^2)+mast+(1|focal) +(1|dyad)+ (1|year),family="binomial")
summary(m1)
plot(ggpredict(m1,terms = c("logDist","relationshipFine")),facet=TRUE,add.data=TRUE)
plot(simulateResiduals(m1, plot = F)) # basic validation check
#printOutput(m1,"Nov_maleSurvival")

# full grid model
m2 <- glmmTMB(data=fullGrid[sex=="M" & distance>0],survival~relationshipFine*poly(logDist,2)+age+I(age^2)+mast+(1|focal)+(1|dyad)+ (1|year),family="binomial")
summary(m2)
plot(simulateResiduals(m2, plot = F)) # basic validation check
#printOutput(m2,"Nov_maleSurvival_fullGrid")



# Male reproductive success -----------------------------------------------

# color for visualization
vcols <- mcols
# calculate sample sizes
fgNo0 <- fullGrid[distance>0 & !is.na(maleARS),,]
fgNo0[,numByRelation:=.N,by=c("sex","relationshipFine")]
fgNo0[,group_col:=relationshipFine,]
sampleSizes <- fgNo0[fgNo0[ , .I[sample(.N,1)] , by = c("sex","group_col")]$V1]

# non-linear visualization for full grid dataset
m <- gam(data=fullGrid[sex=="M" & distance>0],maleARS ~ relationshipFine + s(logDist,by=relationshipFine,k=3)+age+I(age^2)+mast,poisson)
summary.gam(m)
#png("./figures/Nov_final_FigX_maleARSGAMoutput.png",width=9,height=7,units='in',res=700)
plotGAM_viridis(m,sampleSizes,"M","Influence of nearby kin on male reproductive success","Pups sired",20,2.7,20)
#dev.off()

# social neighborhood model
m1 <- glmmTMB(data=fullNH[sex=="M" & distance>0],maleARS~logDist*relationshipFine+age+I(age^2)+mast+(1|focal) +(1|dyad)+ (1|year),family="poisson")
summary(m1)
plot(simulateResiduals(m1, plot = F)) # validation, marginally significant overdispersion but not major
DHARMa::testDispersion(m1)
# zero-inflated version
m1z <- glmmTMB(data=fullNH[sex=="M" & distance>0],maleARS~logDist*relationshipFine+age+I(age^2)+mast+(1|focal) +(1|dyad)+ (1|year),ziformula = ~1,family="poisson")
summary(m1z)
plot(simulateResiduals(m1z, plot = F)) # validation
#printOutput(m1z,"Nov_maleARS")


# full grid model
m2 <- glmmTMB(data=fullGrid[sex=="M" & distance>0],maleARS~relationshipFine*poly(logDist,2)+age+I(age^2)+mast+(1|focal)+(1|dyad)+ (1|year),family="poisson")
summary(m2)
plot(simulateResiduals(m2, plot = F)) # validation
DHARMa::testDispersion(m2) # non-significant overdispersion test
plot(simulateResiduals(m2, plot = F)) # validation
#printOutput(m2,"Nov_maleARS_fullGrid")












