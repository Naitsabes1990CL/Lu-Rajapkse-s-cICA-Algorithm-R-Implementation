###Connect to localhost###
library("RMySQL", lib.loc="~/R/win-library/3.5")
#Connect to data base
con<-dbConnect(RMySQL::MySQL(), host = "localhost",dbname="electrical_distribution_ds",user = "root", password = "")

#Fetch fragments
table_f1=dbSendQuery(con, "SELECT A.Date, A.Hour, A.ID, A.Measurement FROM f1_raw_data_2013 as A")
f1=fetch(table_f1, n=-1)
table_f2=dbSendQuery(con, "SELECT A.Date, A.Hour, A.ID, A.Measurement FROM f2_raw_data_2013 as A")
f2=fetch(table_f2, n=-1)
table_f3=dbSendQuery(con, "SELECT A.Date, A.Hour, A.ID, A.Measurement FROM f3_raw_data_2013 as A")
f3=fetch(table_f3, n=-1)
table_f4=dbSendQuery(con, "SELECT A.Date, A.Hour, A.ID, A.Measurement FROM f4_raw_data_2013 as A")
f4=fetch(table_f4, n=-1)
table_f5=dbSendQuery(con, "SELECT A.Date, A.Hour, A.ID, A.Measurement FROM f5_raw_data_2013 as A")
f5=fetch(table_f5, n=-1)

df=rbind.data.frame(as.data.frame(f1), as.data.frame(f2), as.data.frame(f3), as.data.frame(f4), as.data.frame(f5))

###Data Cleaning and visualizations###
library("lubridate", lib.loc="~/R/win-library/3.5")
library("tidyverse", lib.loc="~/R/win-library/3.5")
df['Date']<-as.Date(as.character(df$Date), format=c("%d/%m/%Y"))
df$DateTime<-paste0(as.character(df$Date)," ",as.character(df$Hour),":00:00")
df['DateTime']<-as.POSIXct(df$DateTime, format=c("%Y-%m-%d %H:%M:%S"), tz="GMT")

#Exploring extreme values and outliers
library("tsoutliers", lib.loc="~/R/win-library/3.5")
library("forecast", lib.loc="~/R/win-library/3.5")
new_df<-df %>% select(ID, DateTime, Measurement) 

ID_SET<-unique(df$ID)
df_clean<-data.frame()
df_final<-data.frame()
for (i in ID_SET){
aux<-filter(new_df, ID==i)
measurement_clean<-tsclean(as.ts(aux$Measurement), replace.missing=TRUE, lambda="auto")
measurement_clean<-abs(measurement_clean)
df_clean<-data.frame(aux$ID, aux$DateTime, as.numeric(measurement_clean))
df_final<-rbind.data.frame(df_final,df_clean)
  }

colnames(df_final)<-c('ID', 'DateTime', 'Measurement')
detach("package:tidyverse", unload=TRUE)
detach("package:tidyr", unload=TRUE)
library("reshape2", lib.loc="~/R/win-library/3.5")
SubEstDataSet=dcast(df_final, DateTime  ~ ID, fun.aggregate=mean, value.var = "Measurement")
setwd("C:/Users/sorel/Desktop/Paper ICA/Scripts")
write.csv(SubEstDataSet, file = "SubEstDataSet.csv", sep=",")



