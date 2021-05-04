library(dplyr)
library(keras)

data_12<-inner_join(otu_normalized[1:992,],unique(all_city_samples_table[,2:3]),by=c('city'='group'))
data_12_1<-inner_join(data_12,data_all,by=c('group2'='city'))
row.names(data_12_1)<-colnames(data_used_12)
#### sampling test ####
data_training<-data_12_1
data_test<-otu_normalized[993:1108,]

deep_multi_mystery<-function(data_training,data_test,n,iter,k){
  time_start<-Sys.time()
  ## detail table 
  output_summer_temp_detail<-as.data.frame(matrix(0,nrow=nrow(data_test),ncol=iter+1))
  output_winter_temp_detail<-as.data.frame(matrix(0,nrow=nrow(data_test),ncol=iter+1))
  output_summer_humidity_detail<-as.data.frame(matrix(0,nrow=nrow(data_test),ncol=iter+1))
  output_winter_humidity_detail<-as.data.frame(matrix(0,nrow=nrow(data_test),ncol=iter+1))
  output_population_score_12_detail<-as.data.frame(matrix(0,nrow=nrow(data_test),ncol=iter+1))
  output_coastal_detail<-as.data.frame(matrix(0,nrow=nrow(data_test),ncol=iter+1))
  ## colnames
  colnames(output_summer_temp_detail)<-c(paste('iter',1:iter,sep='_'),'vote')
  colnames(output_winter_temp_detail)<-c(paste('iter',1:iter,sep='_'),'vote')
  colnames(output_summer_humidity_detail)<-c(paste('iter',1:iter,sep='_'),'vote')
  colnames(output_winter_humidity_detail)<-c(paste('iter',1:iter,sep='_'),'vote')
  colnames(output_population_score_12_detail)<-c(paste('iter',1:iter,sep='_'),'vote')
  colnames(output_coastal_detail)<-c(paste('iter',1:iter,sep='_'),'vote')
  ## rownames
  row.names(output_summer_temp_detail)<-row.names(data_test)
  row.names(output_winter_temp_detail)<-row.names(data_test)
  row.names(output_summer_humidity_detail)<-row.names(data_test)
  row.names(output_winter_humidity_detail)<-row.names(data_test)
  row.names(output_population_score_12_detail)<-row.names(data_test)
  row.names(output_coastal_detail)<-row.names(data_test)
  
  for(i in 1:iter){
    set.seed(i)
    ind<-rep(0,nrow(data_training))
    ind[sampling(data_training,k)]<-1
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
  
  result_for_test_detail<-list(output_summer_temp_detail=output_summer_temp_detail,
                               output_winter_temp_detail=output_winter_temp_detail,
                               output_summer_humidity_detail=output_summer_humidity_detail,
                               output_winter_humidity_detail=output_winter_humidity_detail,
                               output_population_score_12_detail=output_population_score_12_detail,
                               output_coastal_detail=output_coastal_detail)
  
  return(result_for_test_detail)
  time_end<-Sys.time()
  time_diff<-time_end-time_start
  print(time_diff)
}

get_mode_number<-function(data_input){
  f_table<-table(as.matrix(data_input))
  f_table_order<-f_table[order(f_table,decreasing = T)]
  return(f_table_order[1])
}
get_second_vote<-function(data_input){
  f_table<-table(as.matrix(data_input))
  f_table_order<-f_table[order(f_table,decreasing = T)]
  return(names(f_table_order)[2])
}
get_second_vote_number<-function(data_input){
  f_table<-table(as.matrix(data_input))
  f_table_order<-f_table[order(f_table,decreasing = T)]
  return(f_table_order[2])
}

weather_output_mystery$output_summer_temp_detail$vote_number<-apply(weather_output_mystery$output_summer_temp_detail[,1:100],1,get_mode_number)
weather_output_mystery$output_summer_temp_detail$vote_2<-apply(weather_output_mystery$output_summer_temp_detail[,1:100],1,get_second_vote)
weather_output_mystery$output_summer_temp_detail$vote_2_number<-apply(weather_output_mystery$output_summer_temp_detail[,1:100],1,get_second_vote_number)
weather_output_mystery$output_summer_temp_detail$vote_2[1]<-0

weather_output_mystery$output_winter_temp_detail$vote_number<-apply(weather_output_mystery$output_winter_temp_detail[,1:100],1,get_mode_number)
weather_output_mystery$output_winter_temp_detail$vote_2<-apply(weather_output_mystery$output_winter_temp_detail[,1:100],1,get_second_vote)
weather_output_mystery$output_winter_temp_detail$vote_2_number<-apply(weather_output_mystery$output_winter_temp_detail[,1:100],1,get_second_vote_number)

weather_output_mystery$output_summer_humidity_detail$vote_number<-apply(weather_output_mystery$output_summer_humidity_detail[,1:100],1,get_mode_number)
weather_output_mystery$output_summer_humidity_detail$vote_2<-apply(weather_output_mystery$output_summer_humidity_detail[,1:100],1,get_second_vote)
weather_output_mystery$output_summer_humidity_detail$vote_2_number<-apply(weather_output_mystery$output_summer_humidity_detail[,1:100],1,get_second_vote_number)
weather_output_mystery$output_summer_humidity_detail$vote_2[33]<-1

weather_output_mystery$output_winter_humidity_detail$vote_number<-apply(weather_output_mystery$output_winter_humidity_detail[,1:100],1,get_mode_number)
weather_output_mystery$output_winter_humidity_detail$vote_2<-apply(weather_output_mystery$output_winter_humidity_detail[,1:100],1,get_second_vote)
weather_output_mystery$output_winter_humidity_detail$vote_2_number<-apply(weather_output_mystery$output_winter_humidity_detail[,1:100],1,get_second_vote_number)

weather_output_mystery$output_population_score_12_detail$vote_number<-apply(weather_output_mystery$output_population_score_12_detail[,1:100],1,get_mode_number)
weather_output_mystery$output_population_score_12_detail$vote_2<-apply(weather_output_mystery$output_population_score_12_detail[,1:100],1,get_second_vote)
weather_output_mystery$output_population_score_12_detail$vote_2_number<-apply(weather_output_mystery$output_population_score_12_detail[,1:100],1,get_second_vote_number)

weather_output_mystery$output_coastal_detail$vote_number<-apply(weather_output_mystery$output_coastal_detail[,1:100],1,get_mode_number)
weather_output_mystery$output_coastal_detail$vote_2<-apply(weather_output_mystery$output_coastal_detail[,1:100],1,get_second_vote)
weather_output_mystery$output_coastal_detail$vote_2_number<-apply(weather_output_mystery$output_coastal_detail[,1:100],1,get_second_vote_number)
weather_output_mystery$output_coastal_detail$vote_2[78]<-0
















