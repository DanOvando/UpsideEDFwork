#install.packages("rfishbase", repos = c("http://packages.ropensci.org", "http://cran.rstudio.com"), type="source")


library(rfishbase)
list_fields("TLinf")

setwd("/Volumes/COMPATIBLE/ContractWork/EDF/Upside/Models/FishBase/")

SpeciesListFB<-load_taxa()
#write.csv(SpeciesListFB, "SpeciesListFB.csv", row.names=F)
SpecListFB<-paste(as.character(SpeciesListFB$Genus), as.character(SpeciesListFB$Species))
SpecListFB2<-as.vector(SpecListFB)

# check<-list_fields("Temp")
# head(check)
# check2<-data.frame(check)
# TL_data<-ecology(SpecListFB2[1:200], fields=c("DietTLu", "DietseTLu"))
# 
# home_data<-species(SpecListFB2[1:200])
# write.csv(home_data, "home_data_FB.csv", row.names=F)
setwd("/Volumes/COMPATIBLE/ContractWork/EDF/Upside/Models/FishBase/FullFishBase_April2016/")

mat<-maturity(SpecListFB2)
mat<-data.frame(mat)
MaturityData<-mat[,c(2:3,6:7,9)]
write.csv(mat, "maturityFB.csv", row.names=F)
write.csv(MaturityData, "MaturityData.csv", row.names=F)
avetm<-aggregate()

popgwth<-popgrowth(SpecListFB2)
popgwth<-data.frame(popgwth)
names(popgwth)
TempGrowth<-popgwth[,c(2:3,9,17:22,67:68)]
write.csv(popgwth, "popgrowthFB.csv", row.names=F)
write.csv(TempGrowth, "TempGrowthData.csv", row.names=F)

species_Length<-species(SpecListFB2)
species_Length<-data.frame(species_Length)
names(species_Length)
LengthData<-species_Length[,c(1,37,47)]
write.csv(species_Length, "species.csv", row.names=F)
write.csv(LengthData, "LengthData.csv", row.names=F)

stocks_Temp<-stocks(SpecListFB2)
stocks_Temp<-data.frame(stocks_Temp)
names(stocks_Temp)
StocksTempData<-stocks_Temp[,c(2,25:26, 28)]
write.csv(stocks_Temp, "stock.csv", row.names=F)
write.csv(StocksTempData, "StocksTempData.csv", row.names=F)

#####################################
setwd('/Volumes/COMPATIBLE/ContractWork/EDF/Upside/Models/FishBase/FullFishBase_April2016')
MaturityData<-read.csv("MaturityData.csv", stringsAsFactors=F, header=T)
head(MaturityData)
MaturityData$tm<-ifelse(is.na(MaturityData$tm), MaturityData$AgeMatMin, MaturityData$tm)
MaturityData$tm<-ifelse(is.na(MaturityData$tm), MaturityData$AgeMatMin2, MaturityData$tm)
avetm<-aggregate(tm~sciname, data=MaturityData, mean)
write.csv(avetm, "avetm.csv", row.names=T)

LengthData<-read.csv("LengthData.csv", stringsAsFactors=F, header=T)

StocksTempData<-read.csv("StocksTempData.csv", stringsAsFactors=F, header=T)
StocksTempData$Temperature2<-rowMeans(StocksTempData[,c("TempMin", "TempMax")], na.rm=TRUE)
StocksTempData$Temperature2<-ifelse(StocksTempData$Temperature2=="NaN", NA, StocksTempData$Temperature2)
StocksTempData2<-aggregate(Temperature2~sciname, data=StocksTempData, mean)

TempGrowthData<-read.csv("TempGrowthData.csv", stringsAsFactors=F, header=T)
TempGrowthData$TLinfinity<-ifelse(is.na(TempGrowthData$TLinfinity), TempGrowthData$Loo, TempGrowthData$TLinfinity)
TempGrowthData$phi1<-2*log10(TempGrowthData$TLinfinity)+log10(TempGrowthData$K)
names(TempGrowthData)
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

SpeciesListFB$sciname<-paste(SpeciesListFB$Genus, SpeciesListFB$Species)
Species_mpack<-merge(SpeciesListFB, avetm, by="sciname", all.x=T)
Species_mpack<-merge(Species_mpack, LengthData, by="sciname", all.x=T)
Species_mpack<-merge(Species_mpack, aveVar, by="sciname", all.x=T)
Species_mpack<-merge(Species_mpack, StocksTempData2, by="sciname", all.x=T)
Species_mpack$Temperature3<-ifelse(is.na(Species_mpack$Temperature)==T, Species_mpack$Temperature2, Species_mpack$Temperature)
write.csv(Species_mpack, "Species_mpack.csv", row.names=F)



#####################################

lw<-length_weight(SpecListFB2[1:200])
species("Sardinops sagax", field="Length")
mat<-maturity("Sardinops sagax", field="tm")
matave<-mean(mat$tm, na.rm=T)
#popgwth<-popgrowth("Sardinops sagax")
write.csv(popgwth, "popgrwth_SM.csv", row.names=F)
pg_M<-mean(popgwth$K) # this gets you what's in mpack
TLinfK<-popgrowth("Sardinops sagax", field=c("TLinfinity", "K"))
write.csv(TLinfK, "TLinfK.csv", row.names=F)
#prob dont't want the following:
#plf<-poplf("Sardinella brasiliensis", field="GrowthK")
temp1<-popgrowth("Sardinella brasiliensis", field="Temperature")
temp1_m<-mean(temp1$Temperature, na.rm=T)
tempstock<-stocks("Sardinella brasiliensis", field=c("TempMin", "TempMax"))
tempstock_m<-mean(tempstock$TempMin, tempstock$TempMax, na.rm=T)
#temp2<-popqb("Sardinops sagax", field="Temperature")
#temp2_m<-mean(temp2$Temperature, na.rm=T)



lw2<-data.frame(lw)
names(lw2)
write.csv(age, "Pomacanthus maculosus_LW.csv", row.names=F)
write.csv(lw, "length_weight_FB.csv", row.names=F)

growth<-popgrowth(SpecListFB2[1:200])
write.csv(growth, "growth_FB.csv", row.names=F)

age<-popgrowth(SpecListFB2[1:200])
age<-popgrowth("Pomacanthus maculosus")

write.csv(age, "Pomacanthus maculosus.csv", row.names=F)

options(FISHBASE_API = "http://fishbase.ropensci.org")
SpeciesListSLB<-load_taxa(server="http://fishbase.ropensci.org/sealifebase")
write.csv(SpeciesListSLB, "SpeciesListSLB.csv", row.names=F)
SpeciesListSLB<-paste(as.character(SpeciesListSLB$Genus), as.character(SpeciesListSLB$Species))
SpeciesListSLB2<-as.vector(SpeciesListSLB)

check<-stocks("Cetorhinus maximus", c("TempMin", "TempMax", "StockDefs"))
check<-data.frame(check)
write.csv(check, "check.csv", row.names=F)


