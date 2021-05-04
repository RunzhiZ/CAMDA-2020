## Bogota
table_bogota<-mystery_table[which(mystery_table$true=='bogota'),]
table_bogota$check<-0
table_bogota[which(table_bogota$pred_1=='mystery'|table_bogota$pred_2=='mystery'),'check']<-1

## Hong Kong
table_HKG<-mystery_table[which(mystery_table$true=='HKG'),]
table_HKG$check<-0
table_HKG[which(table_HKG$pred_1=='HKG'|table_HKG$pred_2=='HKG'),'check']<-1

## IEV 
table_IEV<-mystery_table[which(mystery_table$true=='IEV'),]
table_IEV$check<-0
table_IEV[which(table_IEV$pred_1=='IEV'|table_IEV$pred_2=='IEV'),'check']<-1

## karkow
table_krakow<-mystery_table[which(mystery_table$true=='krakow'),]
table_krakow$check<-0
table_krakow[which(table_krakow$pred_1=='mystery'|table_krakow$pred_2=='mystery'),'check']<-1

## marseille
table_marseille<-mystery_table[which(mystery_table$true=='marseille'),]
table_marseille$check<-0
table_marseille[which(table_marseille$pred_1=='mystery'|table_marseille$pred_2=='mystery'),'check']<-1

## naples
table_naples<-mystery_table[which(mystery_table$true=='naples'),]
table_naples$check<-0
table_naples[which(table_naples$pred_1=='mystery'|table_naples$pred_2=='mystery'),'check']<-1

## TPE
table_TPE<-mystery_table[which(mystery_table$true=='TPE'),]
table_TPE$check<-0
table_TPE[which(table_TPE$pred_1=='TPE'|table_TPE$pred_2=='TPE'),'check']<-1

## TYO
table_TYO<-mystery_table[which(mystery_table$true=='TYO'),]
table_TYO$check<-0
table_TYO[which(table_TYO$pred_1=='TYO_17'|table_TYO$pred_2=='TYO_17'),'check']<-1

## Vienna
table_vienna<-mystery_table[which(mystery_table$true=='vienna'),]
table_vienna$check<-0
table_vienna[which(table_vienna$pred_1=='mystery'|table_vienna$pred_2=='mystery'),'check']<-1

## ZRH
table_ZRH<-mystery_table[which(mystery_table$true=='ZRH'),]
table_ZRH$check<-0
table_ZRH[which(table_ZRH$pred_1=='ZRH'|table_ZRH$pred_2=='ZRH'),'check']<-1

1-mean(rbind(table_bogota$check,
           table_krakow$check,
           table_marseille$check,
           table_naples$check,
           table_vienna$check))

#### weather ####
mystery_weather<-cbind(row.names(weather_output_mystery$output_summer_temp_detail),
                       weather_output_mystery$output_summer_temp_detail[,101],
                       weather_output_mystery$output_winter_temp_detail[,101],
                       weather_output_mystery$output_summer_humidity_detail[,101],
                       weather_output_mystery$output_winter_humidity_detail[,101],
                       weather_output_mystery$output_population_score_12_detail[,101],
                       weather_output_mystery$output_coastal_detail[,101]
                       )
colnames(mystery_weather)<-c('sample','summer_temp','winter_temp','summer_humidity','winter_humidity',
                             'population_score_12','coastal')
mystery_weather<-as.data.frame(mystery_weather)
mystery_weather[,2:7]<-apply(mystery_weather[,2:7],2,as.numeric)

table_bogota<-inner_join(table_bogota,mystery_weather,by=c('sample'='sample'))
table_HKG<-inner_join(table_HKG,mystery_weather,by=c('sample'='sample'))
table_IEV<-inner_join(table_IEV,mystery_weather,by=c('sample'='sample'))
table_krakow<-inner_join(table_krakow,mystery_weather,by=c('sample'='sample'))
table_marseille<-inner_join(table_marseille,mystery_weather,by=c('sample'='sample'))
table_naples<-inner_join(table_naples,mystery_weather,by=c('sample'='sample'))
table_TPE<-inner_join(table_TPE,mystery_weather,by=c('sample'='sample'))
table_TYO<-inner_join(table_TYO,mystery_weather,by=c('sample'='sample'))
table_vienna<-inner_join(table_vienna,mystery_weather,by=c('sample'='sample'))
table_ZRH<-inner_join(table_ZRH,mystery_weather,by=c('sample'='sample'))



colMeans(table_bogota[which(table_bogota$pred_1=='mystery'|table_bogota$pred_2=='mystery'),10:15])
colMeans(table_bogota[,10:15])
colMeans(table_HKG[,10:15])
colMeans(table_IEV[,10:15])
colMeans(table_krakow[,10:15])
colMeans(table_krakow[which(table_krakow$pred_1=='mystery'|table_krakow$pred_2=='mystery'),10:15])
colMeans(table_marseille[,10:15])
colMeans(table_naples[,10:15])
colMeans(table_naples[which(table_naples$pred_1=='mystery'|table_naples$pred_2=='mystery'),10:15])

colMeans(table_TPE[,10:15])
colMeans(table_TYO[,10:15])
colMeans(table_vienna[,10:15])
colMeans(table_vienna[which(table_vienna$pred_1=='mystery'|table_vienna$pred_2=='mystery'),10:15])

colMeans(table_ZRH[,10:15])

length(which(mystery_table$true=='bogota'|mystery_table$true=='krakow'|mystery_table$true=='vienna'|mystery_table$true=='naples'|mystery_table$true=='marseille'))
length(which((mystery_table$true=='bogota'|mystery_table$true=='krakow'|mystery_table$true=='vienna'|mystery_table$true=='naples'|mystery_table$true=='marseille')&(mystery_table$pred_1=='mystery'|mystery_table$pred_2=='mystery')))
