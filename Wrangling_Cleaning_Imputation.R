###DATA CLEANING/WRANGLING###
library("lubridate", lib.loc="~/R/win-library/3.5")
library("tidyverse", lib.loc="~/R/win-library/3.5")
library("padr", lib.loc="~/R/win-library/3.5")
library("forecast", lib.loc="~/R/win-library/3.5")
setwd("C:/Users/sorel/Desktop/Paper ICA/Data Chilectra/2013/SNAPS 2013/SNAPS Avg Power")

#1. Load Feeder Data Set
AlimDataSet=read.csv("MasterAlims2013DataSet.csv", header = TRUE, sep = ",")

#2. Generate Valid Time-stamp
AlimDataSet['Fecha']<-as.Date(as.character(AlimDataSet$Fecha), format=c("%d/%m/%Y"))
AlimDataSet$DateTime<-paste0(as.character(AlimDataSet$Fecha)," ",as.character(AlimDataSet$H),":00:00")
AlimDataSet['DateTime']<-as.POSIXct(AlimDataSet$DateTime, format=c("%Y-%m-%d %H:%M:%S"), tz="GMT")

#Sanity Check 1: How many unique days? (ie: less than 365 -> Timestamp imputation is required
#length(unique(as.Date(AlimDataSet$DateTime, format=c("%d/%m/%Y"))))

#3. Select relevant columns
MasterAlimDataSet=data.frame(AlimDataSet%>%dplyr::select(ID, AlimID, DateTime, PTotal)%>%group_by(ID, AlimID, DateTime))
names(MasterAlimDataSet)<-c('SSEE_ID', 'Alim_ID', 'DateTime', 'Measurement')


#Sanity Check 2: Check for Duplicates (ie: they were removed in the data ingestion python script)
#DuplicatesAlim=MasterAlimDataSet[duplicated(MasterAlimDataSet),]
#head(DuplicatesAlim)

#4. Explore NA's
AlimNA=MasterAlimDataSet%>%filter(is.na(Measurement)==TRUE)
#Which Feeders and to what substation have NA's measurements?
data.frame(unique(AlimNA$Alim_ID), unique(AlimNA$SSEE_ID))
Alim125=MasterAlimDataSet%>%filter(Alim_ID==125)
Alim341=MasterAlimDataSet%>%filter(Alim_ID==341)
Alim319=MasterAlimDataSet%>%filter(Alim_ID==319)

#5. Explore Date Completeness: SHOULD be unique 365 days and 367 Feeders ID
length(unique(as.Date(MasterAlimDataSet$DateTime, format=c("%d/%m/%Y"))))
length(unique(MasterAlimDataSet$Alim_ID)) #Check Length of Alim_IDs

#6.From LONG to WIDE format
detach("package:tidyverse", unload=TRUE)
detach("package:tidyr", unload=TRUE)
library("reshape2", lib.loc="~/R/win-library/3.5")
M=dcast(MasterAlimDataSet, DateTime  ~ Alim_ID, fun.aggregate=mean, value.var = "Measurement")
detach("package:reshape2", unload=TRUE)

#7. Add Missing Dates
library("padr", lib.loc="~/R/win-library/3.5")
library("tidyverse", lib.loc="~/R/win-library/3.5")
FeedersID=paste0("Alim_",unique(MasterAlimDataSet$Alim_ID))
MasterAlimDF=data.frame(pad(M, interval="hour"))
names(MasterAlimDF)=c('DateTime', FeedersID)

#8. Impute Missing Values
library("Amelia", lib.loc="~/R/win-library/3.5")
library("mice", lib.loc="~/R/win-library/3.5")
library("naniar", lib.loc="~/R/win-library/3.5")
setwd("C:/Users/sorel/Desktop/Paper ICA/Scripts")

#9. Compute Missing Value Percentage
ColNA=data.frame(t(data.frame(MasterAlimDF[,2:length(MasterAlimDF)] %>% dplyr::select(everything()) 
                              %>% summarise_all(funs(sum(is.na(.))))))*(1/dim(MasterAlimDF)[1])*100)
ColNA=cbind(rownames(ColNA), data.frame(ColNA, row.names=NULL))
names(ColNA)=c('Alim_ID','Missing_Values')

#10. Generate Visualizations Distribution AND Concentration
pMissing1<-vis_miss(MasterAlimDF[,2:length(MasterAlimDF)], warn_large_data = FALSE)+ 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave(pMissing1, file="./Missing Value Plots/Missing Values Distribution.PNG", width = 14, height = 10, units = "cm")
pMissing2<-gg_miss_var(MasterAlimDF, show_pct = TRUE)+
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
ggsave(pMissing2, file="./Missing Value Plots/Missing Values Concentration.PNG", width = 14, height = 10, units = "cm")

#11. Closer Inspection to Missing Values Feeders
#a. More than 10% of the observations are missing values per feeder
HighMissingValues=ColNA%>%filter(Missing_Values>10)
MissingValuesDataSet=MasterAlimDF%>%dplyr::select(DateTime, Alim_54, Alim_111,  
                                           Alim_122, Alim_125, Alim_126, Alim_165, Alim_233, Alim_259, Alim_341, Alim_347)

#Sanity Check 3: Plot them and analyze their functional form.  
#HighMVAlimsID=unique(as.character(HighMissingValues$Alim_ID))
#for (i in 1:length(HighMVAlimsID)){
#  pMVAlim<-ggplot(MissingValuesDataSet, aes(x = DateTime, y = MissingValuesDataSet[,i+1])) + geom_miss_point() + theme_dark()
#  file_name = paste("./Missing Value Plots/MVAlim_", HighMVAlimsID[i], ".PNG", sep="")
#  ggsave(pMVAlim, file=file_name, width = 14, height = 10, units = "cm")}

#Sanity Check 4:Contrast with Sample of Low Missing Values Feeders
#LowMissingValues=ColNA%>%filter(Missing_Values<=10)
#LowMVAlimsID=unique(as.character(LowMissingValues$Alim_ID))
#random_index=sample(1:357, 10)
#SampleIDs=LowMVAlimsID[random_index]
#SampleLowMVDataSet=MasterAlimDF%>%dplyr::select(DateTime, starts_with(SampleIDs[1]))
#for (i in 2:length(random_index)){
#  TempLowMVDataSet=MasterAlimDF%>%dplyr::select(DateTime, starts_with(SampleIDs[i]))
#  SampleLowMVDataSet=data.frame(SampleLowMVDataSet,TempLowMVDataSet[,2])
#}
#names(SampleLowMVDataSet)=c('DateTime', SampleIDs)
#
#for (i in 1:length(SampleIDs)){
#  pLowMVAlim<-ggplot(SampleLowMVDataSet, aes(x = DateTime, y = SampleLowMVDataSet[,i+1])) + geom_miss_point() + theme_dark()
#  file_name = paste("./Missing Value Plots/LowMVAlim_", SampleIDs[i], ".PNG", sep="")
#  ggsave(pLowMVAlim, file=file_name, width = 14, height = 10, units = "cm")
#  }

#b. Less than 10% of observations are missing values per feeder
LowMissingValues=ColNA%>%filter(Missing_Values<=10)
LowMVAlimsID=unique(as.character(LowMissingValues$Alim_ID))
LowMVDataSet=MasterAlimDF%>%dplyr::select(DateTime, starts_with(LowMVAlimsID[1]))
for (i in 2:length(LowMVAlimsID)){
  TempLowMVDataSet=MasterAlimDF%>%dplyr::select(DateTime, starts_with(LowMVAlimsID[i]))
  LowMVDataSet=data.frame(LowMVDataSet,TempLowMVDataSet[,2])
}
names(LowMVDataSet)=c('DateTime', LowMVAlimsID)

#c. Check Distributions of missing values for months, weekdays and hours.
DatesNA=data.frame(LowMVDataSet[rowSums(is.na(LowMVDataSet)) > 0,]%>%dplyr::select(DateTime)%>%group_by(DateTime))
unique(months(DatesNA[,1])) #Which months
unique(weekdays(DatesNA[,1])) #Which weekdays
sort(as.numeric(unique(strftime(DatesNA[,1], format="%H"))), decreasing=FALSE)

#13. Impute Missing Values using Mixed Chained Equations with Predictive Mean Matching (Non-Parametric Approach)
library(mice)
tempMiceImputed=mice(data.frame(LowMVDataSet[,2:dim(LowMVDataSet)[2]]))
MiceImputed=data.frame(LowMVDataSet$DateTime,complete(tempMiceImputed))
names(MiceImputed)=c("DateTime", colnames(LowMVDataSet[,2:dim(LowMVDataSet)[2]]))
write.csv(MiceImputed, file = "MiceImputed.csv", sep=",")









