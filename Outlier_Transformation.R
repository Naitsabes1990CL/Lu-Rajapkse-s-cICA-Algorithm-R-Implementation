###DATA TRANSFORMATION AND OUTLIER DETECTION###

#1. Load MICE Imputed Data
library("dplyr", lib.loc="~/R/win-library/3.5")
setwd("C:/Users/sorel/Desktop/Paper ICA/Scripts")
MiceImputed=data.frame(read.csv("MiceImputed.csv", header = TRUE, sep = ","))
MiceImputed=MiceImputed[,3:dim(MiceImputed)[2]]
DataFeeders=MiceImputed[,2:dim(MiceImputed)[2]]

##########################################################2. OUTLIER FLAGGING PIPELINE##########################################################
#2.1 LOAD FACTOR OUTLIERS
library("lubridate", lib.loc="~/R/win-library/3.5")
library("reshape2", lib.loc="~/R/win-library/3.5")
library("ggplot2", lib.loc="~/R/win-library/3.5")
library("Amelia", lib.loc="~/R/win-library/3.5")
library("mice", lib.loc="~/R/win-library/3.5")
library("naniar", lib.loc="~/R/win-library/3.5")
DataFeeders.LF=DataFeeders
DataFeeders.LF$Date=as.Date(MiceImputed[,1])
DataFeeders.LF$Week=week(as.POSIXct(MiceImputed[,1]))
##############################################################2.1.1 DAILY LOAD FACTOR##########################################################
#Long Format#
DataFeedersLong=melt(DataFeeders.LF,
                     # ID variables - all the variables to keep but not split apart on
                     id.vars=c("Date", "Week"),
                     # The source columns
                     measure.vars=colnames(DataFeeders.LF[,1:355]),
                     # Name of the destination column that will identify the original
                     # column that the measurement came from
                     variable.name="Alim_ID",
                     value.name="Measurement"
)
DailyMeanDataFeeders=data.frame(DataFeedersLong %>% select(Date, Alim_ID, Measurement)
                                %>%group_by(Date, Alim_ID) %>% summarise(Mean = mean(abs(Measurement), na.rm=TRUE)))
DailyMaxDataFeeders=data.frame(DataFeedersLong %>% select(Date, Alim_ID, Measurement)
                               %>% group_by(Date, Alim_ID) %>% summarise(Max = max(abs(Measurement), na.rm=TRUE)))
DailyLoadFactor=data.frame(DailyMeanDataFeeders$Date, DailyMeanDataFeeders$Alim_ID, (DailyMeanDataFeeders$Mean/DailyMaxDataFeeders$Max))
names(DailyLoadFactor)=c("Date", "Alim_ID", "Load.Factor")

#Wide Format#
LoadFactorDataFeeders=dcast(DailyLoadFactor, Date  ~ Alim_ID, fun.aggregate=mean, value.var = "Load.Factor")
PctNA=data.frame(t(data.frame(LoadFactorDataFeeders[,2:length(LoadFactorDataFeeders)] %>% dplyr::select(everything()) 
                              %>% summarise_all(funs(sum(is.na(.))))))*(1/dim(LoadFactorDataFeeders)[1])*100)
PctNA$Alim_ID=row.names(PctNA)
names(PctNA)=c("pct.NA", "Alim_ID")
Daily.LF.Outliers=PctNA%>%dplyr::filter(pct.NA>0)%>%dplyr::select(Alim_ID)

#2.1.2 WEEKLY LOAD FACTOR
#Long Format#
WeeklyMeanDataFeeders=data.frame(DataFeedersLong %>% select(Week, Alim_ID, Measurement)
                                 %>%group_by(Week, Alim_ID) %>% summarise(Mean = mean(abs(Measurement), na.rm=TRUE)))
WeeklyMaxDataFeeders=data.frame(DataFeedersLong %>% select(Week, Alim_ID, Measurement)
                                %>% group_by(Week, Alim_ID) %>% summarise(Max = max(abs(Measurement), na.rm=TRUE)))
WeeklyLoadFactor=data.frame(WeeklyMeanDataFeeders$Week, WeeklyMeanDataFeeders$Alim_ID, (WeeklyMeanDataFeeders$Mean/WeeklyMaxDataFeeders$Max))
names(WeeklyLoadFactor)=c("Week", "Alim_ID", "Load.Factor")
#Wide Format#
LoadFactorDataFeeders2=dcast(WeeklyLoadFactor, Week  ~ Alim_ID, fun.aggregate=mean, value.var = "Load.Factor")
PctNA2=data.frame(t(data.frame(LoadFactorDataFeeders2[,2:length(LoadFactorDataFeeders2)] %>% dplyr::select(everything()) 
                               %>% summarise_all(funs(sum(is.na(.))))))*(1/dim(LoadFactorDataFeeders)[1])*100)
PctNA2$Alim_ID=row.names(PctNA2)
names(PctNA2)=c("pct.NA2", "Alim_ID")
Weekly.LF.Outliers=PctNA2%>%dplyr::filter(pct.NA2>0)%>%dplyr::select(Alim_ID)

######2.1.3 OUTPUT OF LOAD FACTOR OUTLIERS######
Load.Factor.Outliers=as.character(t(left_join(Daily.LF.Outliers, Weekly.LF.Outliers, by=c("Alim_ID", "Alim_ID"))))
idx=match(Load.Factor.Outliers, names(DataFeeders))
Output.LF.Outliers=DataFeeders[,-idx]

######################################################2.2 FREQUENCY/PERIOD OUTLIERS###############################################################
library("aTSA", lib.loc="~/R/win-library/3.5")
library("TSA", lib.loc="~/R/win-library/3.5")


#2.2.1 Generate Periodgrams using spec.gram(unbiased with taper=0.1)
p=spec.pgram(Output.LF.Outliers[,1], plot=FALSE, taper=0.1, demean=TRUE)
freq=p$freq[which(p$spec>=sort(p$spec, decreasing = TRUE)[2], arr.ind = TRUE)]#Find top 3 values
M=rbind.data.frame(freq)
names(M)=c("Top1", "Top2")
#Iterate over the set of feeders
for (i in 2:dim(Output.LF.Outliers)[2]){
  p=spec.pgram(Output.LF.Outliers[,i], plot=FALSE, taper=0.1, demean=TRUE)
  freq=p$freq[which(p$spec>=sort(p$spec, decreasing = TRUE)[2], arr.ind = TRUE)]#Find top 3 values
  M=rbind.data.frame(M, freq)
}
M=data.frame(colnames(Output.LF.Outliers), M)
names(M)=c("Alim_ID", "Top1", "Top2")
M$ID=as.numeric(gsub("[a-zA-Z_ ]", "", M$Alim_ID))

#2.2.2 Transform M into daily periods 1/f=Period
PM=data.frame(M[,1], M[,4], round(1/M[,2]), round(1/M[,3]))
names(PM)=c("Alim_ID", "ID", "Period.1", "Period.2")
#Hourly in a year: 1/(365*24)=0.0001141553
#Hourly in a weeek: 1/(7*24)=0.005952381

#2.2.3 Compute Period Anomalies
Period.1.counts=data.frame(PM %>% count(Period.1, sort = TRUE))
Period.2.counts=data.frame(PM %>% count(Period.2, sort = TRUE))

#2.2.4 Flag Period Outliers
PM.1=left_join(PM, Period.1.counts, by = c("Period.1" = "Period.1"))
names(PM.1)=c("Alim_ID", "ID", "Period.1", "Period.2", "n.Period.1")
PM.2=left_join(PM.1, Period.2.counts, by=c("Period.2" = "Period.2"))
names(PM.2)=c("Alim_ID", "ID", "Period.1", "Period.2", "n.Period.1", "n.Period.2")
Period.outliers=PM.2%>%dplyr::filter(n.Period.1<(dim(Output.LF.Outliers)[2]*0.1) & n.Period.2<(dim(Output.LF.Outliers)[2]*0.1))%>%select(Alim_ID)

######2.2.5 OUTPUT OF PERIOD OUTLIERS######
Anomally.Period.Outliers=as.character(t(Period.outliers))
idx=match(Anomally.Period.Outliers, names(Output.LF.Outliers))
Output.Period.Outliers=Output.LF.Outliers[,-idx]

######################################################################2.3 UN-BALANCED CLUSTER OUTLIERS############################################
library("dplyr", lib.loc="~/R/win-library/3.5")
library("tsoutliers", lib.loc="~/R/win-library/3.5")
library("TSclust", lib.loc="~/R/win-library/3.5")
library("factoextra", lib.loc="~/R/win-library/3.5")
library("dendextend", lib.loc="~/R/win-library/3.5")
library("cluster", lib.loc="C:/Program Files/R/R-3.5.1/library")
library("NbClust", lib.loc="~/R/win-library/3.5")

#2.3.1 Compute Dissimilarity Matrix
TSMatrix2=t(as.matrix(Output.Period.Outliers))
DissMat2=diss(TSMatrix2, "INT.PER")

#2.3.2 Select Number of Clusters
#DoParallel routine
require("doParallel")
# Create parallel workers
workers <- makeCluster(6L)
# Preload dtwclust in each worker; not necessary but useful
invisible(clusterEvalQ(workers, library("factoextra")))
# Register the backend; this step MUST be done
registerDoParallel(workers)
set.seed(12345)
AvgSilhouette=fviz_nbclust(Output.Period.Outliers, FUN = hcut, hc_func=c("hclust"), hc_method=c("complete"), method = "silhouette")
# Stop parallel workers
stopCluster(workers)
# Go back to sequential computation
registerDoSEQ()
#Extract Maximum Silhouette
Cluster.Sil=AvgSilhouette$data
CLUS.NUM=as.numeric(Cluster.Sil[which.max(Cluster.Sil$y),1])

#2.3.2 Hierarchical Aglomerative Clustering Complete Linkage
hcPER.complete=hclust(DissMat2, method="complete")
plot(hcPER.complete)
rect.hclust(hcPER.complete, k=CLUS.NUM, border=2:12) 
sub_grp.complete=cutree(hcPER.complete, k = CLUS.NUM)
t=data.frame(sub_grp.complete)
Output.HClust=data.frame(row.names(t), t, row.names=NULL)
names(Output.HClust)=c("Alim_ID", "Cluster")
#2.3.3 Construct Summary DataFrame
ClusterSum=data.frame(Output.HClust%>%dplyr::count(Cluster))
#2.3.4 Filter by less than 10% of observations
Outlier.Clust=ClusterSum%>%dplyr::filter(n<dim(DissMat2)[2]*0.1)%>%dplyr::select(Cluster)
#2.3.5 Flag Alim_ID associated only with Outlier.Clust
H.Clust.Outliers=data.frame()
for (i in 1:dim(Outlier.Clust)[1]){
H.Clust.Outliers=rbind.data.frame(H.Clust.Outliers,Output.HClust%>%filter(Cluster==Outlier.Clust[i,1])%>%select(Alim_ID))
}
######2.3.5 OUTPUT OF UN-BALANCED CLUSTER OUTLIERS######
Un.Balanced.Cluster.Outliers=as.character(t(H.Clust.Outliers))
idx=match(Un.Balanced.Cluster.Outliers, names(Output.Period.Outliers))
Output.Un.Balanced.Cluster.Outliers=Output.Period.Outliers[,-idx]
CleanedDataSet=data.frame(MiceImputed[,1],Output.Un.Balanced.Cluster.Outliers)
write.csv(CleanedDataSet, file = "CleanedFeeder.csv", sep=",")

###################################################3. Transform Series into MultipleSeasonal Time Series#########################################
#Multiple Seasonal Adjustment for hourly and weekly seasonalities for TS with frequency=3600 (hourly measurements)
library("aTSA", lib.loc="~/R/win-library/3.5")
library("tseries", lib.loc="~/R/win-library/3.5")
library("TSA", lib.loc="~/R/win-library/3.5")
library("forecast", lib.loc="~/R/win-library/3.5")

#3.1 Noise (Season Adjustment substracting the seasonal pattern)
SDS=data.frame(seasadj(stl(msts(Output.Un.Balanced.Cluster.Outliers[,1], seasonal.periods=c(24,7*24)), "periodic")))
for (i in 2:dim(Output.Un.Balanced.Cluster.Outliers)[2]){
  SDS=data.frame(SDS, seasadj(stl(msts(Output.Un.Balanced.Cluster.Outliers[,i], seasonal.periods=c(24,7*24)), "periodic")))
}
SDS=data.frame(MiceImputed[,1], SDS)
names(SDS)=c("DateTime", paste0("Noise_",colnames(Output.Un.Balanced.Cluster.Outliers)))
write.csv(SDS, file = "Noise_Term.csv", sep=",")

#3.2 Seasonal Pattern (Substracting the noice from the obtained TS)
SeasonD=data.frame(Output.Un.Balanced.Cluster.Outliers[,1]-seasadj(stl(msts(Output.Un.Balanced.Cluster.Outliers[,1], 
                                                                            seasonal.periods=c(24,7*24)), "periodic")))
for (i in 2:dim(Output.Un.Balanced.Cluster.Outliers)[2]){
  SeasonD=data.frame(SeasonD, Output.Un.Balanced.Cluster.Outliers[,i]-seasadj(stl(msts(Output.Un.Balanced.Cluster.Outliers[,i], 
                                                                                       seasonal.periods=c(24,7*24)), "periodic")))
}
SeasonD=data.frame(MiceImputed[,1], SeasonD)
names(SeasonD)=c("DateTime", paste0("Season_",colnames(Output.Un.Balanced.Cluster.Outliers)))
write.csv(SeasonD, file = "Seasonal_Term.csv", sep=",")

#3.3 Noise+Seasonal Data Set
Noise_Season_df=cbind.data.frame(SeasonD, SDS[,-1])
write.csv(Noise_Season_df, file = "Noise_Season_Output.csv", sep=",")




