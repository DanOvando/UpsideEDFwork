######################################
#--------------------------------------------------
#This code pulls current fishbase life history data from the server and assembles a dataframe
##Created by: Kristin Kleisner
##Date: April 14, 2016

###TO DO: need to assemble this into a list like mpack???
######################################

###Load latest copy of rfishbase from ropensci and call library:
#install.packages("rfishbase", repos = c("http://packages.ropensci.org", "http://cran.rstudio.com"), type="source")
library(rfishbase)

###Set working directory:
setwd("/Volumes/COMPATIBLE/ContractWork/EDF/Upside/Models/FishBase/Upside_FishBase_git")

###Use the following command to search through all FishBase tables for particular fields---sometimes multiple tables will have data:
list_fields("TLinf")


###Load the full list of taxa in FishBase (>33000 species):
SpeciesListFB<-load_taxa()
#write.csv(SpeciesListFB, "SpeciesListFB.csv", row.names=F)
SpecListFB<-paste(as.character(SpeciesListFB$Genus), as.character(SpeciesListFB$Species))
SpecListFB2<-as.vector(SpecListFB)

###The following code (commented out) pulls the entire table of data for all available species base on the full FishBase species list. The full datasets are saved as well as the subset of columns containing data for the mpack assembly.
###NOTE: These pulls from FishBase take a couple of hours potentially, so I can send these files to be read in below!!!

# mat<-maturity(SpecListFB2)
# mat<-data.frame(mat)
# MaturityData<-mat[,c(2:3,6:7,9)]
# write.csv(mat, "maturityFB.csv", row.names=F)
# write.csv(MaturityData, "MaturityData.csv", row.names=F)
# avetm<-aggregate()
# 
# popgwth<-popgrowth(SpecListFB2)
# popgwth<-data.frame(popgwth)
# names(popgwth)
# TempGrowth<-popgwth[,c(2:3,9,17:22,67:68)]
# write.csv(popgwth, "popgrowthFB.csv", row.names=F)
# write.csv(TempGrowth, "TempGrowthData.csv", row.names=F)
# 
# species_Length<-species(SpecListFB2)
# species_Length<-data.frame(species_Length)
# names(species_Length)
# LengthData<-species_Length[,c(1,37,47)]
# write.csv(species_Length, "species.csv", row.names=F)
# write.csv(LengthData, "LengthData.csv", row.names=F)
# 
# stocks_Temp<-stocks(SpecListFB2)
# stocks_Temp<-data.frame(stocks_Temp)
# names(stocks_Temp)
# StocksTempData<-stocks_Temp[,c(2,25:26, 28)]
# write.csv(stocks_Temp, "stock.csv", row.names=F)
# write.csv(StocksTempData, "StocksTempData.csv", row.names=F)

#####################################
###Set working directory where the datatables are located and read in files:
setwd('/Volumes/COMPATIBLE/ContractWork/EDF/Upside/Models/FishBase/Upside_FishBase_git/Data')

###Read in dataset with maturity data for species:
MaturityData<-read.csv("MaturityData.csv", stringsAsFactors=F, header=T)
head(MaturityData)
MaturityData$tm<-ifelse(is.na(MaturityData$tm), MaturityData$AgeMatMin, MaturityData$tm)
MaturityData$tm<-ifelse(is.na(MaturityData$tm), MaturityData$AgeMatMin2, MaturityData$tm)
avetm<-aggregate(tm~sciname, data=MaturityData, mean)

###Read in dataset with length data for species:
LengthData<-read.csv("LengthData.csv", stringsAsFactors=F, header=T)

###Read in dataset with temperature data for stocks:
StocksTempData<-read.csv("StocksTempData.csv", stringsAsFactors=F, header=T)
StocksTempData$Temperature2<-rowMeans(StocksTempData[,c("TempMin", "TempMax")], na.rm=TRUE)
StocksTempData$Temperature2<-ifelse(StocksTempData$Temperature2=="NaN", NA, StocksTempData$Temperature2)
StocksTempData2<-aggregate(Temperature2~sciname, data=StocksTempData, mean)

###Read in dataset with temperature data and growth model parameters:
TempGrowthData<-read.csv("TempGrowthData.csv", stringsAsFactors=F, header=T)
TempGrowthData$TLinfinity<-ifelse(is.na(TempGrowthData$TLinfinity), TempGrowthData$Loo, TempGrowthData$TLinfinity)
TempGrowthData$phi1<-2*log10(TempGrowthData$TLinfinity)+log10(TempGrowthData$K)

###The following script calculates the average K from all models for a species based on Loo and phi-prime:
aveTemp<-aggregate(Temperature~sciname, data=TempGrowthData, mean)
aveK<-aggregate(K~sciname, data=TempGrowthData, mean)
names(aveK)<-c("sciname", "aveK")
aveTLinf<-aggregate(TLinfinity~sciname, data=TempGrowthData, mean)
aveTLinf$logTlinf<-2*log10(aveTLinf$TLinfinity)
names(aveTLinf)<-c("sciname", "aveTLinf", "logTLinf")
avephi1<-aggregate(phi1~sciname, data=TempGrowthData, mean)
names(avephi1)<-c("sciname", "avephi1")
aveVar<-merge(aveK, aveTLinf, by="sciname", all.x=T)
aveVar<-merge(aveVar, avephi1, by="sciname", all.x=T)
aveVar$aveKphi<-10^(aveVar$avephi1-aveVar$logTLinf)
aveVar<-merge(aveVar, aveTemp, by="sciname", all.x=T)
#write.csv(aveVar, "aveVar.csv", row.names=F)


###Merge datasets into one dataframe for all species--this preserves the family info in case we want to calculate family level values to fill in for  species with missing life history info (i.e., fill in with family level data):
SpeciesListFB$sciname<-paste(SpeciesListFB$Genus, SpeciesListFB$Species)
Species_mpack<-merge(SpeciesListFB, avetm, by="sciname", all.x=T)
Species_mpack<-merge(Species_mpack, LengthData, by="sciname", all.x=T)
Species_mpack<-merge(Species_mpack, aveVar, by="sciname", all.x=T)
Species_mpack<-merge(Species_mpack, StocksTempData2, by="sciname", all.x=T)
Species_mpack$Temperature3<-ifelse(is.na(Species_mpack$Temperature)==T, Species_mpack$Temperature2, Species_mpack$Temperature)
#write.csv(Species_mpack, "Species_mpack.csv", row.names=F)



