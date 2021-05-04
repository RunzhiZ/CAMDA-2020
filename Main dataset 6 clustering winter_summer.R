library(dplyr)
library(cluster)
library(tidyverse)
library(factoextra)
library(NbClust)

load('Main dataset 5 city data processing')
#### Summer ####
summercols<-grep("Summer",colnames(data_weather_season_2),value=TRUE)
summerweather<-data_weather_season_2%>%select(all_of(summercols))
rownames(summerweather)<-data_weather_season_2$City

## Summer Temp ##
summer_temp<-grep("tempC",colnames(summerweather),value=TRUE)
summer_tempweather<-data_weather_season_2%>%select(all_of(summer_temp))
summer_tempweather<-summer_tempweather[,1:6]
NbClust(summer_tempweather,method="kmeans")
# fviz_nbclust(summer_tempweather,kmeans,method="silhouette")
# fviz_nbclust(summer_tempweather,kmeans,method="wss")
# gap_stat<-clusGap(summer_tempweather,FUN=kmeans,nstart=25,K.max=10,B=100)
# fviz_gap_stat(gap_stat)
k_summer_temp<-kmeans(summer_tempweather,centers=4,nstart=25)
fviz_cluster(k_summer_temp, data = summer_tempweather,main='Summer temp')

## Summer Humidity ##
summer_humidity<-grep("humidity",colnames(summerweather),value=TRUE)
summer_humidityweather<-data_weather_season_2%>%select(all_of(summer_humidity))

NbClust(summer_humidityweather,method="kmeans")
# fviz_nbclust(summer_humidityweather,kmeans,method="silhouette")
# fviz_nbclust(summer_humidityweather,kmeans,method="wss")
# gap_stat<-clusGap(summer_humidityweather,FUN=kmeans,nstart=25,K.max=10,B=100)
# fviz_gap_stat(gap_stat)
k_summer_humidity<-kmeans(summer_humidityweather,centers=3,nstart=25)
fviz_cluster(k_summer_humidity, data = summer_humidityweather,main='Summer humidity')

#### winter ####
wintercols<-grep("Winter",colnames(data_weather_season_2),value=TRUE)
winterweather<-data_weather_season_2%>%select(all_of(wintercols))
rownames(summerweather)<-data_weather_season_2$City

## winter Temp ##
winter_temp<-grep("tempC",colnames(winterweather),value=TRUE)
winter_tempweather<-data_weather_season_2%>%select(all_of(winter_temp))
winter_tempweather[19,8]<-23.5
NbClust(winter_tempweather,method="kmeans")
# fviz_nbclust(winter_tempweather,kmeans,method="silhouette")
# fviz_nbclust(winter_tempweather,kmeans,method="wss")
# gap_stat<-clusGap(winter_tempweather,FUN=kmeans,nstart=25,K.max=10,B=100)
# fviz_gap_stat(gap_stat)
k_winter_temp<-kmeans(winter_tempweather,centers=5,nstart=25)
fviz_cluster(k_winter_temp, data = winter_tempweather,main='winter_temp')

## winter Humidity ##
winter_humidity<-grep("humidity",colnames(winterweather),value=TRUE)
winter_humidityweather<-data_weather_season_2%>%select(all_of(winter_humidity))
NbClust(winter_humidityweather,method="kmeans")
fviz_nbclust(winter_humidityweather,kmeans,method="silhouette")
fviz_nbclust(winter_humidityweather,kmeans,method="wss")
# gap_stat<-clusGap(winter_humidityweather,FUN=kmeans,nstart=25,K.max=10,B=100)
# fviz_gap_stat(gap_stat)
k_winter_humidity<-kmeans(winter_humidityweather,centers=2,nstart=25)
fviz_cluster(k_winter_humidity, data = winter_humidityweather,main='Winter_humidity')

#### categorical data frame ####
data_all<-as.data.frame(matrix(NA,nrow=23,ncol=7))
colnames(data_all)<-c('city','summer_temp','winter_temp','summer_humidity','winter_humidity','population_score','coastal')
data_all$city<-data_weather_season_2$City_Code
data_all$summer_temp<-k_summer_temp$cluster
data_all$winter_temp<-k_winter_temp$cluster
data_all$summer_humidity<-k_summer_humidity$cluster
data_all$winter_humidity<-k_winter_humidity$cluster
data_all$population_score<-data_biomes$population_score
data_all$coastal<-data_location$Coastal_City

data_all[which(k_summer_temp$cluster==3),'summer_temp']<-2
data_all[which(k_summer_temp$cluster==4),'summer_temp']<-2
data_all[which(k_summer_temp$cluster==1),'summer_temp']<-0
data_all[which(k_summer_temp$cluster==2),'summer_temp']<-1

data_all[which(k_summer_humidity$cluster==3),'summer_humidity']<-1
data_all[which(k_summer_humidity$cluster==1),'summer_humidity']<-2
data_all[which(k_summer_humidity$cluster==2),'summer_humidity']<-0

data_all[which(k_winter_temp$cluster==3),'winter_temp']<-3
data_all[which(k_winter_temp$cluster==4),'winter_temp']<-2
data_all[which(k_winter_temp$cluster==2),'winter_temp']<-1
data_all[which(k_winter_temp$cluster==1),'winter_temp']<-0
data_all[which(k_winter_temp$cluster==5),'winter_temp']<-0

data_all[which(k_winter_humidity$cluster==1),'winter_humidity']<-0
data_all[which(k_winter_humidity$cluster==2),'winter_humidity']<-1

data_all[data_all$population_score<=11,'population_score']<-0
data_all[data_all$population_score>=12,'population_score']<-1
data_all$coastal<-as.numeric(data_all$coastal)
data_all$coastal<-as.numeric(data_all$coastal)

save.image('Main dataset 6 clustering winter_summer.RData')








