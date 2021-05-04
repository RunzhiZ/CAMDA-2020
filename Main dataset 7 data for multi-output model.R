library(dplyr)
library(keras)
setwd('E:/University of Florida/My research/CAMDA 2020/data')
colnames(data_all)[6]<-'population_score_12'

data_12<-inner_join(data_used_normalized_12,unique(all_city_samples_table[,2:3]),by=c('city'='group'))
data_12_1<-inner_join(data_12,data_all,by=c('group2'='city'))
row.names(data_12_1)<-row.names(data_used_normalized_12)
#### sampling test ####
set.seed(123)
test_ind<-sampling(data_used_normalized_12,5)
data_12_training<-data_12_1[-test_ind,]
data_12_test<-data_12_1[test_ind,]

#### multi-output ####
deep_multi_f<-function(data_training,data_test,n,iter,k,data_training_sampling){
  time_start<-Sys.time()
  ## detail table 
  output_summer_temp_detail<-as.data.frame(matrix(0,nrow=nrow(data_test),ncol=iter+3))
  output_winter_temp_detail<-as.data.frame(matrix(0,nrow=nrow(data_test),ncol=iter+3))
  output_summer_humidity_detail<-as.data.frame(matrix(0,nrow=nrow(data_test),ncol=iter+3))
  output_winter_humidity_detail<-as.data.frame(matrix(0,nrow=nrow(data_test),ncol=iter+3))
  output_population_score_12_detail<-as.data.frame(matrix(0,nrow=nrow(data_test),ncol=iter+3))
  output_coastal_detail<-as.data.frame(matrix(0,nrow=nrow(data_test),ncol=iter+3))
  ## colnames
  colnames(output_summer_temp_detail)<-c(paste('iter',1:iter,sep='_'),'weather_data','vote','check')
  colnames(output_winter_temp_detail)<-c(paste('iter',1:iter,sep='_'),'weather_data','vote','check')
  colnames(output_summer_humidity_detail)<-c(paste('iter',1:iter,sep='_'),'weather_data','vote','check')
  colnames(output_winter_humidity_detail)<-c(paste('iter',1:iter,sep='_'),'weather_data','vote','check')
  colnames(output_population_score_12_detail)<-c(paste('iter',1:iter,sep='_'),'weather_data','vote','check')
  colnames(output_coastal_detail)<-c(paste('iter',1:iter,sep='_'),'weather_data','vote','check')
  ## rownames
  row.names(output_summer_temp_detail)<-row.names(data_test)
  row.names(output_winter_temp_detail)<-row.names(data_test)
  row.names(output_summer_humidity_detail)<-row.names(data_test)
  row.names(output_winter_humidity_detail)<-row.names(data_test)
  row.names(output_population_score_12_detail)<-row.names(data_test)
  row.names(output_coastal_detail)<-row.names(data_test)
  
  output_summer_temp_detail[,'weather_data']<-data_test[,'summer_temp']
  output_winter_temp_detail[,'weather_data']<-data_test[,'winter_temp']
  output_summer_humidity_detail[,'weather_data']<-data_test[,'summer_humidity']
  output_winter_humidity_detail[,'weather_data']<-data_test[,'winter_humidity']
  output_population_score_12_detail[,'weather_data']<-data_test[,'population_score_12']
  output_coastal_detail[,'weather_data']<-data_test[,'coastal']
  
  for(i in 1:iter){
    set.seed(i)
    ind<-rep(0,nrow(data_training))
    ind[sampling(data_training_sampling,k)]<-1
    data_deep_training <- data_training[ind==0,]

    ## input
    main_input<-layer_input(shape=c(n),name='main_input')
    second_layer_input<-main_input %>%
      layer_dense(units = 128, activation = 'relu') %>% 
      layer_dense(units = 128, activation = 'relu') %>%
      layer_dropout(rate=0.5) %>%
      layer_dense(units = 128, activation = 'relu') %>%
      layer_dropout(rate=0.5)
    ## output
    output_summer_temp<-second_layer_input %>%
      layer_dense(units=3,activation = 'softmax',name='output_summer_temp')
    output_winter_temp<-second_layer_input %>%
      layer_dense(units=4,activation = 'softmax',name='output_winter_temp')
    output_summer_humidity<-second_layer_input %>%
      layer_dense(units=3,activation = 'softmax',name='output_summer_humidity')
    output_winter_humidity<-second_layer_input %>%
      layer_dense(units=1,activation = 'sigmoid',name='output_winter_humidity')
    output_population_score_12<-second_layer_input %>%
      layer_dense(units=1,activation = 'sigmoid',name='output_population_score_12')
    output_coastal<-second_layer_input %>%
      layer_dense(units=1,activation = 'sigmoid',name='output_coastal')
    
    model<-keras_model(
      inputs=main_input,
      outputs = c(output_summer_temp,
                  output_winter_temp,
                  output_summer_humidity,
                  output_winter_humidity,
                  output_population_score_12,
                  output_coastal)
    )
    ## model
    model %>% compile(
      optimizer = 'rmsprop',
      loss = list(output_summer_temp = 'categorical_crossentropy',
                  output_winter_temp = 'categorical_crossentropy',
                  output_summer_humidity = 'categorical_crossentropy',
                  output_winter_humidity = 'binary_crossentropy',
                  output_population_score_12 = 'binary_crossentropy',
                  output_coastal = 'binary_crossentropy'),
      metrics = 'accuracy'
    )
    model %>% fit(
      x = list(main_input = as.matrix(data_deep_training[,1:n])),
      y = list(output_summer_temp = to_categorical(data_deep_training[,47]), 
               output_winter_temp = to_categorical(data_deep_training[,48]),
               output_summer_humidity = to_categorical(data_deep_training[,49]), 
               output_winter_humidity = data_deep_training[,50],
               output_population_score_12 = data_deep_training[,51],
               output_coastal = data_deep_training[,52]),
      epochs = 500,
      verbose=0,
      batch_size = 32
    )
    pred <- model %>% predict(as.matrix(data_test[,1:n]))
    output_summer_temp_detail[,i]<-apply(pred[[1]],1,which.max)-1
    output_winter_temp_detail[,i]<-apply(pred[[2]],1,which.max)-1
    output_summer_humidity_detail[,i]<-apply(pred[[3]],1,which.max)-1
    output_winter_humidity_detail[,i]<-round(pred[[4]])
    output_population_score_12_detail[,i]<-round(pred[[5]])
    output_coastal_detail[,i]<-round(pred[[6]])
  }
  output_summer_temp_detail$vote<-apply(output_summer_temp_detail[,1:iter],1,getmode)
  output_winter_temp_detail$vote<-apply(output_winter_temp_detail[,1:iter],1,getmode)
  output_summer_humidity_detail$vote<-apply(output_summer_humidity_detail[,1:iter],1,getmode)
  output_winter_humidity_detail$vote<-apply(output_winter_humidity_detail[,1:iter],1,getmode)
  output_population_score_12_detail$vote<-apply(output_population_score_12_detail[,1:iter],1,getmode)
  output_coastal_detail$vote<-apply(output_coastal_detail[,1:iter],1,getmode)
  
  for(i in 1:nrow(output_summer_temp_detail)){
    if(output_summer_temp_detail[i,'vote']==output_summer_temp_detail[i,'weather_data']){output_summer_temp_detail[i,'check']<-1}
  }
  for(i in 1:nrow(output_winter_temp_detail)){
    if(output_winter_temp_detail[i,'vote']==output_winter_temp_detail[i,'weather_data']){output_winter_temp_detail[i,'check']<-1}
  }
  for(i in 1:nrow(output_summer_humidity_detail)){
    if(output_summer_humidity_detail[i,'vote']==output_summer_humidity_detail[i,'weather_data']){output_summer_humidity_detail[i,'check']<-1}
  }
  for(i in 1:nrow(output_winter_humidity_detail)){
    if(output_winter_humidity_detail[i,'vote']==output_winter_humidity_detail[i,'weather_data']){output_winter_humidity_detail[i,'check']<-1}
  }
  for(i in 1:nrow(output_population_score_12_detail)){
    if(output_population_score_12_detail[i,'vote']==output_population_score_12_detail[i,'weather_data']){output_population_score_12_detail[i,'check']<-1}
  }
  for(i in 1:nrow(output_coastal_detail)){
    if(output_coastal_detail[i,'vote']==output_coastal_detail[i,'weather_data']){output_coastal_detail[i,'check']<-1}
  }
  
  output_summer_temp_table<-output_summer_temp_detail[,c('weather_data','vote','check')]
  output_winter_temp_table<-output_winter_temp_detail[,c('weather_data','vote','check')]
  output_summer_humidity_table<-output_summer_humidity_detail[,c('weather_data','vote','check')]
  output_winter_humidity_table<-output_winter_humidity_detail[,c('weather_data','vote','check')]
  output_population_score_12_table<-output_population_score_12_detail[,c('weather_data','vote','check')]
  output_coastal_table<-output_coastal_detail[,c('weather_data','vote','check')]
  
  result_for_test_detail<-list(output_summer_temp_detail=output_summer_temp_detail,
                               output_winter_temp_detail=output_winter_temp_detail,
                               output_summer_humidity_detail=output_summer_humidity_detail,
                               output_winter_humidity_detail=output_winter_humidity_detail,
                               output_population_score_12_detail=output_population_score_12_detail,
                               output_coastal_detail=output_coastal_detail)
  result_for_test<-list(output_summer_temp_table=output_summer_temp_table,
                        output_winter_temp_table=output_winter_temp_table,
                        output_summer_humidity_table=output_summer_humidity_table,
                        output_winter_humidity_table=output_winter_humidity_table,
                        output_population_score_12_table=output_population_score_12_table,
                        output_coastal_table=output_coastal_table)
  
  time_end<-Sys.time()
  time_diff<-time_end-time_start
  print(time_diff)
  result_output<-list(result_for_test_detail=result_for_test_detail,
                      result_for_test=result_for_test,
                      time=time_diff)
  return(result_output)
  
}


#### multi output ####
weather_output_12<-deep_multi_f(data_12_training,data_12_test,44,100,5,data_training_12)

#### weather for test ####
data_12_test_weather<-as.data.frame(matrix(NA,nrow=199))
data_12_test_weather$summer_temp<-weather_output_12_3$result_for_test[[1]]$vote
data_12_test_weather$winter_temp<-weather_output_12_3$result_for_test[[2]]$vote
data_12_test_weather$summer_humidity<-weather_output_12_3$result_for_test[[3]]$vote
data_12_test_weather$winter_humidity<-weather_output_12_3$result_for_test[[4]]$vote
data_12_test_weather$population_score_12<-weather_output_12_3$result_for_test[[5]]$vote
data_12_test_weather$coastal<-weather_output_12_3$result_for_test[[6]]$vote
data_12_test_weather<-data_12_test_weather[2:7]


data_12_test_weather<-cbind(data_12_test_weather[,c(-45,-46)],data_12_test_weather[,45])
colnames(data_12_test_weather)[51]<-'city'
data_12_training_weather<-cbind(data_12_training[,c(-45,-46)],data_12_training[,45])
colnames(data_12_training_weather)[51]<-'city'























