library(dplyr)
#### read data ####
data_biomes<-read.csv("biomes_meta.csv",header = T,check.names = F)
data_weather<-read.csv("weather_meta.csv",header = T,check.names = F)
data_location<-read.csv("locations_meta.csv",header = T,check.names = F)
#### use only winter/summer data ####
data_weather_season<-select(data_weather,contains(c('City','December','January','February','June','July','August')))
data_weather_season_1<-select(data_weather_season,contains(c('City','_tempC_avg','_maxtempC_max','_mintempC_min',
                                                             'humidity_avg','humidity_max','humidity_min',
                                                             'pressure_avg','pressure_max','pressure_min')))
data_weather_season_2<-data_weather_season_1
colnames(data_weather_season_2)<-gsub('December','Winter_1',colnames(data_weather_season_2))
colnames(data_weather_season_2)<-gsub('January','Winter_2',colnames(data_weather_season_2))
colnames(data_weather_season_2)<-gsub('February','Winter_3',colnames(data_weather_season_2))
colnames(data_weather_season_2)<-gsub('June','Summer_1',colnames(data_weather_season_2))
colnames(data_weather_season_2)<-gsub('July','Summer_2',colnames(data_weather_season_2))
colnames(data_weather_season_2)<-gsub('August','Summer_3',colnames(data_weather_season_2))

#### transform data for cities from southern 
data_weather_north_summer<-select(data_weather_season_2,contains(c('Summer')))[15:16,]
data_weather_north_winter<-select(data_weather_season_2,contains(c('Winter')))[15:16,]
data_weather_season_2[15:16,colnames(data_weather_north_winter)]<-data_weather_north_summer
data_weather_season_2[15:16,colnames(data_weather_north_summer)]<-data_weather_north_winter

#### average ####
data_weather_season_avg<-as.data.frame(matrix(NA,nrow=23,ncol=14))
colnames(data_weather_season_avg)<-c('City','City_Code','Winter_tempC_avg','Summer_tempC_avg','Winter_tempC_std','Summer_tempC_std',
                                     'Winter_humidity_avg','Summer_humidity_avg','Winter_humidity_std','Summer_humidity_std',
                                     'Winter_pressure_avg','Summer_pressure_avg','Winter_pressure_std','Summer_pressure_std')
data_weather_season_avg[,c('City','City_Code')]<-data_weather_season_2[,1:2]
data_weather_season_avg[,'Winter_tempC_avg']<-apply(data_weather_season_2[,3:5],1,mean)
data_weather_season_avg[,'Summer_tempC_avg']<-apply(data_weather_season_2[,6:8],1,mean)
data_weather_season_avg[,'Winter_tempC_std']<-apply(data_weather_season_2[,9:11],1,mean)
data_weather_season_avg[,'Summer_tempC_std']<-apply(data_weather_season_2[,12:14],1,mean)
data_weather_season_avg[,'Winter_humidity_avg']<-apply(data_weather_season_2[,15:17],1,mean)
data_weather_season_avg[,'Summer_humidity_avg']<-apply(data_weather_season_2[,18:20],1,mean)
data_weather_season_avg[,'Winter_humidity_std']<-apply(data_weather_season_2[,21:23],1,mean)
data_weather_season_avg[,'Summer_humidity_std']<-apply(data_weather_season_2[,24:26],1,mean)
data_weather_season_avg[,'Winter_pressure_avg']<-apply(data_weather_season_2[,27:29],1,mean)
data_weather_season_avg[,'Summer_pressure_avg']<-apply(data_weather_season_2[,30:32],1,mean)
data_weather_season_avg[,'Winter_pressure_std']<-apply(data_weather_season_2[,33:35],1,mean)
data_weather_season_avg[,'Summer_pressure_std']<-apply(data_weather_season_2[,36:38],1,mean)

#### biomes ####
city_biomes<-data.frame(type=c('Dense Settlement Biomes','Village Biomes','Cropland Biomes','Rangeland Biomes','Forested Biomes'),score=5:1)
data_biomes$score_1<-inner_join(data_biomes,city_biomes,by=c('Neighbor_Biome_1_1'='type'))[ncol(data_biomes)+1]
data_biomes$score_2<-inner_join(data_biomes,city_biomes,by=c('Neighbor_Biome_1_2'='type'))[ncol(data_biomes)+1]
data_biomes$score_3<-inner_join(data_biomes,city_biomes,by=c('Neighbor_Biome_1_3'='type'))[ncol(data_biomes)+1]
data_biomes$population_score<-apply(data_biomes[,11:13],1,sum)

#### data_all ####
data_all<-cbind(data_weather_season_avg,data_location[,8],data_biomes[,14])
colnames(data_all)[15:16]<-c('Coastal_city','population_score')

#### save ####
save.image('Main dataset 5 city data processing')






