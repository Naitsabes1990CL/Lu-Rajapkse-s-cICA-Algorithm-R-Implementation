###LOAD SUBEST AND ALIM DATA SET TO CONSTRUCT ADJACENCY MATRIX
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
SSEE.Feeder.DataSet=data.frame(AlimDataSet%>%dplyr::select(ID, AlimID)%>%group_by(ID))
names(SSEE.Feeder.DataSet)=c('SSEE_ID', 'Alim_ID')

#4.Adjacency Matrix: From LONG to WIDE format
detach("package:tidyverse", unload=TRUE)
detach("package:tidyr", unload=TRUE)
library("reshape2", lib.loc="~/R/win-library/3.5")
M=dcast(SSEE.Feeder.DataSet, Alim_ID  ~ SSEE_ID, fun.aggregate=mean, value.var="Alim_ID", fill=-1)
M[,1]=paste0("Alim_",unique(M[,1]))
M.aux=M[,2:dim(M)[2]]
M.aux[M.aux>=0]=1
M.aux[M.aux==-1]=0
Adj.M=data.frame(M[,1], M.aux)
colnames(Adj.M)=c("Alim_ID", paste0("SSEE_",colnames(M.aux)))

###HIERARCHICAL SIMULATION###
#1. Load Both Feeders and SubEst Data Sets
setwd("C:/Users/sorel/Desktop/Paper ICA/Scripts")
#Feeders
Feeders=read.csv("CleanedFeeder.csv", header = TRUE, sep = ",")
Feeders=Feeders[,2:dim(Feeders)[2]]
#Substations
SSEE=read.csv("CleanedSubEst.csv", header = TRUE, sep = ",")
SSEE=SSEE[,2:dim(SSEE)[2]]
colnames(SSEE)=c("DateTime", paste0("SSEE",colnames(SSEE[,2:dim(SSEE)[2]])))
colnames(SSEE)=gsub("X", "_", colnames(SSEE))

#2. Reduce number of Feeders and SubStations to those valid
#SSEE
Valid.SSEE=data.frame(colnames(SSEE[,-1]),1)
colnames(Valid.SSEE)=c("SSEE", "aux")
Adj.M.SSEE=data.frame(colnames(Adj.M[,-1]),1)
#Left Join
colnames(Adj.M.SSEE)=c("SSEE", "aux")
mergeCols=c("SSEE")
left1=merge(Valid.SSEE, Adj.M.SSEE, by = mergeCols, all.x = TRUE, sort=FALSE)
#Output
Final.SSEE=as.character(t(left1%>%dplyr::filter(is.na(aux.y)==FALSE)%>%dplyr::select(SSEE)))

#Feeders
Valid.Feeders=data.frame(colnames(Feeders[,-1]),1)
colnames(Valid.Feeders)=c("Feeders", "aux")
Adj.M.Feeders=data.frame(Adj.M[,1],1)
#Left Join
colnames(Adj.M.Feeders)=c("Feeders", "aux")
mergeCols=c("Feeders")
left2=merge(Valid.Feeders, Adj.M.Feeders, by = mergeCols, all.x = TRUE, sort=FALSE)
#Output
Final.Feeders=as.character(t(left2%>%dplyr::filter(is.na(aux.y)==FALSE)%>%dplyr::select(Feeders)))

#Summation/Adjecency Matrix
col.idx=match(Final.SSEE, names(Adj.M))
row.idx=match(Final.Feeders, Adj.M[,1])
Final.Adj.M=Adj.M[row.idx,col.idx]

#Feeders/SSEE
aux.SSEE=SSEE[,2:dim(SSEE)[2]]
SSEE.Final=aux.SSEE[,match(Final.SSEE, names(aux.SSEE))]
aux.Feeders=Feeders[,2:dim(Feeders)[2]]
Feeders.Final=aux.Feeders[,match(Final.Feeders, names(aux.Feeders))]

#3. Summation Matrix Form
Y=as.matrix(rbind(t(Feeders.Final), t(SSEE.Final)))
A=as.matrix(rbind(t(Final.Adj.M), t(diag(1,325,325))))
S=as.matrix(t(Feeders.Final))

#Extra. ARIMA MODELS FOR FEEDERS
#library("forecast", lib.loc="~/R/win-library/3.5")
#lambda=BoxCox.lambda(as.ts(Feeders[,2]))
#box.cox.feeder.2=BoxCox(as.ts(Feeders[,2]),lambda)
#autoplot(BoxCox(as.ts(Feeders[,2]),lambda))
#fit1=msts(box.cox.feeder.2, seasonal.periods=c(12,24,7*24))
#fit2=auto.arima(as.ts(Feeders[,2]))
#fit3=stl(msts(as.ts(Feeders[,2]), seasonal.periods=c(12,24,7*24)), "periodic") ###WINNER!!!
#forecast::forecast(fit, 365)
#checkresiduals(fit1)
#autoplot(forecast::forecast(fit1, method="rwdrift",365))
#autoplot(forecast::forecast(fit2, 365))
#autoplot(forecast::forecast(fit3$time.series, method="naive", 8760))
#autoplot(forecast::forecast(fit3$time.series[,1], 8760))

#4. Feeders Simulations using SAME mutliseasonal TS model as the data-preprocessing step
Simulated.Feeder=data.frame(Feeders[,1])
for (i in 1:dim(Feeders.Final)[2]){
  fit=stl(msts(as.ts(Feeders.Final[,i]), seasonal.periods=c(12,24,7*24)), "periodic")
  sim=forecast::forecast((fit$time.series[,1]+fit$time.series[,2]+fit$time.series[,3]), h=8760)
  Simulated.Feeder=data.frame(Simulated.Feeder, sim$x)}

#5. Substations Simulations using SAME mutliseasonal TS model as the data-preprocessing step
Simulated.SSEE=data.frame(SSEE[,1])
for (i in 1:dim(SSEE.Final)[2]){
  fit=stl(msts(as.ts(SSEE.Final[,i]), seasonal.periods=c(12,24,7*24)), "periodic")
  sim=forecast::forecast((fit$time.series[,1]+fit$time.series[,2]+fit$time.series[,3]), h=8760)
  Simulated.SSEE=data.frame(Simulated.SSEE, sim$x)}

#6. Save both Simulated Data Sets
write.csv(Simulated.SSEE, file = "Simulated_SSEE.csv", sep=",")
write.csv(Simulated.Feeder, file="Simulated_Feeder.csv", sep=",")



