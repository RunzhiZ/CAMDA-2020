library(dplyr)
library(limma)
library(randomForest)
library(e1071)
library(MASS)
library(ggplot2)
library(edgeR)
library(caret)
library(glmnet)
library(mlbench)
library(keras)

load('Main dataset 4 data for feature selection and machine learning.RData')
colnames(all_city_samples_table)<-c('sample','group','group2')

###### Machine learning results ######
## feature selection
species_0.6<-row.names(all_species_ubiquity_table)[which(all_species_ubiquity_table$ubiquity_0.6>=1)]
family_0.6<-row.names(all_family_ubiquity_table)[which(all_family_ubiquity_table$ubiquity_0.6>=1)]
order_0.6<-row.names(all_order_ubiquity_table)[which(all_order_ubiquity_table$ubiquity_0.6>=1)]
features_0.6_otu_table<-rbind(all_species_table[species_0.6,],all_family_table[family_0.6,],all_order_table[order_0.6,])
#### without weather ####
data_normalized<-normalization_f(features_0.6_otu_table)
#### with weather ####
data_normalized<-inner_join(data_normalized,unique(all_city_samples_table[,2:3]),by=c('city'='group'))
data_normalized<-inner_join(data_normalized,data_all[,-1],by=c('group2'='City_Code'))
data_normalized[,'group2']<-NULL

#### 0.4
result_list_0.4<-list()
for(j in 1:length(city_all)){
  set.seed(123)
  x <- model.matrix(city~., data_normalized)[,-1]
  y <- ifelse(data_normalized$city == city_all[j], 1, 0)
  
  cv.lasso <- cv.glmnet(x, y, alpha = 0.4, family = "binomial")
  # Fit the final model on the training data
  model <- glmnet(x, y, alpha = 0.4, family = "binomial",
                  lambda = cv.lasso$lambda.1se,maxit = 1e+06)
  rm(x,y)
  coef(model)
  result_list_0.4[[j]]<-names(coef(model)[which(coef(model)!=0),])[-1]
}

feature_list_0.4<-as.data.frame(matrix(0,nrow=ncol(data_normalized)-1,ncol=29))
colnames(feature_list_0.4)<-c(city_all,'count')
#### without weather ####
row.names(feature_list_0.4)<-row.names(features_0.6_otu_table)
row.names(feature_list_0.4)[which(row.names(feature_list_0.4)=='JG30-KF-CM45')]<-'JG30.KF.CM45'
#### with weather ####
row.names(feature_list_0.4)<-c(row.names(features_0.6_otu_table),colnames(data_all)[c(-1,-2)])
row.names(feature_list_0.4)[which(row.names(feature_list_0.4)=='JG30-KF-CM45')]<-'JG30.KF.CM45'
row.names(feature_list_0.4)[which(row.names(feature_list_0.4)=='Coastal_city')]<-'Coastal_cityTRUE'

for(i in 1:28){
  feature_list_0.4[result_list_0.4[[i]],i]<-1
}
feature_list_0.4$count<-apply(feature_list_0.4[,1:28],1,sum)

#### Function: sampling for each city ####
sampling<-function(data,n){
  sample_table<-as.data.frame(matrix(NA,nrow=length(table(data$city)),ncol=3))
  colnames(sample_table)<-c('city','number_of_samples','number_for_CV')
  sample_table$city<-names(table(data$city))
  sample_table$number_of_samples<-table(data$city)
  sample_table$number_for_CV<-round(sample_table$number_of_samples/n)
  samples<-NULL
  for(i in 1:length(sample_table$city)){
    samples<-c(samples,sample(which(data$city==sample_table[i,'city']),sample_table[i,'number_for_CV']))
  }
  return(samples)
}

#### Function: get mode ####
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#### Function: get the group ####
tiqu<-function(str){
  substr(str,1,min(which(strsplit(str,split='')[[1]]=='0'))-2)
}
for(i in 1:nrow(all_city_samples_table)){
  all_city_samples_table[i,2]<-tiqu(all_city_samples_table[i,1])
}
colnames(all_city_samples_table)<-c('sample','group','group2')

#### Function: Normalization ####
normalization_f<-function(data){
  data_0<-data[,which(colSums(data)!=0)]
  normalization_t<-data_0
  d0<-DGEList(normalization_t)
  group<-all_city_samples_table[which(all_city_samples_table[,1] %in% colnames(data_0)),2]
  mm<-model.matrix(~0 + group)
  y<-voom(d0, mm, plot = T)$E
  otu_normalized<-as.data.frame(t(y))
  otu_normalized$city<-group
  names(otu_normalized)<-make.names(names(otu_normalized))
  return(otu_normalized)
}
#### Function: deep learning ####
deep_learning_k_fold<-function(data_training,data_test,k1,iter){
  time_start<-Sys.time()
  all_history<-list()
  ## test dataset
  id_table<-data.frame(city=unique(data_training$city),id=0:27)
  data_deep_test<-inner_join(data_test,id_table,by='city')
  data_deep_test<-as.matrix(data_deep_test[,c(1:(ncol(data_deep_test)-2),ncol(data_deep_test))])
  data_deep_test_data<-data_deep_test[, 1:(ncol(data_deep_test)-1)]
  ##data_deep_test_group<-data_deep_test[, 'id']
  ##data_deep_test_Labels<-to_categorical(data_deep_test_group)
  result_for_test_detail<-as.data.frame(matrix(0,nrow=nrow(data_test),ncol=iter+3))
  colnames(result_for_test_detail)<-c(paste('iter',1:iter,sep='_'),'city','vote','check')
  row.names(result_for_test_detail)<-row.names(data_test)
  for(i in 1:nrow(result_for_test_detail)){
    result_for_test_detail[i,'city']<-tiqu(row.names(result_for_test_detail)[i])
  }
  ## training dataset
  for(i in 1:iter){
    set.seed(i)
    ind<-rep(0,nrow(data_training))
    ind[sampling(data_training,k1)]<-1
    data_deep<-inner_join(data_training,id_table,by='city')
    data_deep<-as.matrix(data_deep[,c(1:(ncol(data_deep)-2),ncol(data_deep))])
    data_deep_training <- data_deep[ind==0, 1:(ncol(data_deep)-1)]
    data_deep_validation <- data_deep[ind==1, 1:(ncol(data_deep)-1)]
    data_deep_training_group <- data_deep[ind==0, 'id']
    data_deep_validation_group <- data_deep[ind==1, 'id']
    data_deep_training_Labels <- to_categorical(data_deep_training_group)
    data_deep_validation_Labels <- to_categorical(data_deep_validation_group)
    #### deep learning
    model <- keras_model_sequential() 
    
    # Add layers to the model
    model %>% 
      layer_dense(units = 128, activation = 'relu', input_shape = c(ncol(data_training)-1)) %>% 
      layer_dense(units = 128, activation = 'relu') %>%
      layer_dropout(rate=0.5) %>%
      layer_dense(units = 128, activation = 'relu') %>%
      layer_dropout(rate=0.5) %>%
      layer_dense(units = 28, activation = 'softmax')
    
    # Compile the model
    model %>% compile(
      loss = 'categorical_crossentropy',
      optimizer = 'adam',
      metrics = 'accuracy'
    )
    
    # Store the fitting history in `history` 
    history <- model %>% fit(
      data_deep_training, 
      data_deep_training_Labels, 
      epochs = 300,
      batch_size = 32, 
      verbose=0,
      validation_data=list(data_deep_validation,data_deep_validation_Labels)
    )
    all_history[[i]]<-history
    # Predict the classes for the test data
    classes <- model %>% predict_classes(data_deep_test_data, batch_size = 32)
    # result
    deep_learning_result<-as.data.frame(matrix(0,nrow=length(classes),ncol=2))
    colnames(deep_learning_result)<-c('predict_id','predict')
    deep_learning_result$predict_id<-classes
    deep_learning_result$predict<-inner_join(deep_learning_result,id_table,by=c('predict_id'='id'))[,3]
    result_for_test_detail[,i]<-deep_learning_result$predict
  }
  result_for_test_detail$vote<-apply(result_for_test_detail[,1:iter],1,getmode)
  for(i in 1:nrow(result_for_test_detail)){
    if(result_for_test_detail[i,'vote']==result_for_test_detail[i,'city']){result_for_test_detail[i,'check']<-1}
  }
  result_for_test<-result_for_test_detail[,c('city','vote','check')]
  time_end<-Sys.time()
  time_diff<-time_end-time_start
  print(time_diff)
  result_output<-list(all_history=all_history,
                      result_for_test_detail=result_for_test_detail,
                      result_for_test=result_for_test,
                      time=time_diff)
  return(result_output)
}

#### Function: Random forest - k fold ####
rf_k_fold<-function(data_training,data_test,k1,iter){ 
  time_start<-Sys.time()
  result_for_CV<-list()
  result_for_test_detail<-as.data.frame(matrix(0,nrow=nrow(data_test),ncol=iter+3))
  colnames(result_for_test_detail)<-c(paste('iter',1:iter,sep='_'),'city','vote','check')
  row.names(result_for_test_detail)<-row.names(data_test)
  for(i in 1:nrow(result_for_test_detail)){
    result_for_test_detail[i,'city']<-tiqu(row.names(result_for_test_detail)[i])
  }
  for(i in 1:iter){
    ## Cross-Validation
    set.seed(i)
    k<-sampling(data_training,k1)
    data_train<-data_training[-k,]
    data_validation<-data_training[k,]
    rf_fit<-randomForest(factor(city)~.,data=data_train,ntree=1001,type='classifier')
    rf_pred<-predict(rf_fit,data_validation)
    rf_predict<-as.data.frame(cbind(row.names(data_training)[k],as.vector(rf_pred)))
    colnames(rf_predict)<-c('sample','predict')
    rf_predict$sample<-as.character(rf_predict$sample)
    rf_predict$predict<-as.character(rf_predict$predict)
    rf_predict$city<-0
    for(j in 1:nrow(rf_predict)){
      rf_predict[j,3]<-tiqu(rf_predict[j,1])
    }
    rm(j)
    rf_predict$check<-0
    for(j in 1:nrow(rf_predict)){
      if(rf_predict[j,'predict']==rf_predict[j,'city']){rf_predict[j,'check']<-1}
    }
    rm(j)
    rf_predict<-rf_predict[,c(1,3,2,4)]
    rf_predict_city<-as.data.frame(matrix(0,nrow=length(city_all),ncol=2))
    colnames(rf_predict_city)<-c('city','error_rate')
    for(j in 1:length(city_all)){
      rf_predict_city[j,1]<-city_all[j]
      rf_predict_city[j,2]<-1-mean(rf_predict[which(rf_predict$city==city_all[j]),'check'])
    }
    rm(j)
    rf_predict_city<-rf_predict_city[order(rf_predict_city$error_rate,decreasing = F),]
    error_rate<-1-mean(rf_predict$check)
    result_for_1_iter_CV<-list(error_rate=error_rate,
                               rf_predict=rf_predict,
                               rf_predict_city=rf_predict_city)
    result_for_CV[[i]]<-result_for_1_iter_CV
    
    ## Test prediction
    rf_pred_test<-predict(rf_fit,data_test)
    result_for_test_detail[,i]<-as.vector(rf_pred_test)
  }
  ## vote result
  result_for_test_detail$vote<-apply(result_for_test_detail[,1:iter],1,getmode)
  for(i in 1:nrow(result_for_test_detail)){
    if(result_for_test_detail[i,'vote']==result_for_test_detail[i,'city']){result_for_test_detail[i,'check']<-1}
  }
  result_for_test<-result_for_test_detail[,c('city','vote','check')]
  time_end<-Sys.time()
  time_diff<-time_end-time_start
  print(time_diff)
  result_output<-list(result_for_CV=result_for_CV,
                      result_for_test_detail=result_for_test_detail,
                      result_for_test=result_for_test,
                      time=time_diff)
  return(result_output)
}

#### Function: SVM - k fold ####
svm_k_fold<-function(data_training,data_test,k1,iter){ 
  time_start<-Sys.time()
  result_for_CV<-list()
  result_for_test_detail<-as.data.frame(matrix(0,nrow=nrow(data_test),ncol=iter+3))
  colnames(result_for_test_detail)<-c(paste('iter',1:iter,sep='_'),'city','vote','check')
  row.names(result_for_test_detail)<-row.names(data_test)
  for(i in 1:nrow(result_for_test_detail)){
    result_for_test_detail[i,'city']<-tiqu(row.names(result_for_test_detail)[i])
  }
  for(i in 1:iter){
    ## Cross-Validation
    set.seed(i)
    k<-sampling(data_training,k1)
    data_train<-data_training[-k,]
    data_validation<-data_training[k,]
    svm_fit<-best.svm(factor(city)~.,data=data_train)
    svm_pred<-predict(svm_fit,data_validation)
    svm_predict<-as.data.frame(cbind(row.names(data_training)[k],as.vector(svm_pred)))
    colnames(svm_predict)<-c('sample','predict')
    svm_predict$sample<-as.character(svm_predict$sample)
    svm_predict$predict<-as.character(svm_predict$predict)
    svm_predict$city<-0
    for(j in 1:nrow(svm_predict)){
      svm_predict[j,3]<-tiqu(svm_predict[j,1])
    }
    rm(j)
    svm_predict$check<-0
    for(j in 1:nrow(svm_predict)){
      if(svm_predict[j,'predict']==svm_predict[j,'city']){svm_predict[j,'check']<-1}
    }
    rm(j)
    svm_predict<-svm_predict[,c(1,3,2,4)]
    svm_predict_city<-as.data.frame(matrix(0,nrow=length(city_all),ncol=2))
    colnames(svm_predict_city)<-c('city','error_rate')
    for(j in 1:length(city_all)){
      svm_predict_city[j,1]<-city_all[j]
      svm_predict_city[j,2]<-1-mean(svm_predict[which(svm_predict$city==city_all[j]),'check'])
    }
    rm(j)
    svm_predict_city<-svm_predict_city[order(svm_predict_city$error_rate,decreasing = F),]
    error_rate<-1-mean(svm_predict$check)
    result_for_1_iter_CV<-list(error_rate=error_rate,
                               svm_predict=svm_predict,
                               svm_predict_city=svm_predict_city)
    result_for_CV[[i]]<-result_for_1_iter_CV
    
    ## Test prediction
    svm_pred_test<-predict(svm_fit,data_test)
    result_for_test_detail[,i]<-as.vector(svm_pred_test)
  }
  ## vote result
  result_for_test_detail$vote<-apply(result_for_test_detail[,1:iter],1,getmode)
  for(i in 1:nrow(result_for_test_detail)){
    if(result_for_test_detail[i,'vote']==result_for_test_detail[i,'city']){result_for_test_detail[i,'check']<-1}
  }
  result_for_test<-result_for_test_detail[,c('city','vote','check')]
  time_end<-Sys.time()
  time_diff<-time_end-time_start
  print(time_diff)
  result_output<-list(result_for_CV=result_for_CV,
                      result_for_test_detail=result_for_test_detail,
                      result_for_test=result_for_test,
                      time=time_diff)
  return(result_output)
}

#### Function: LDA - k fold ####
lda_k_fold<-function(data_training,data_test,k1,iter){ 
  time_start<-Sys.time()
  result_for_CV<-list()
  result_for_test_detail<-as.data.frame(matrix(0,nrow=nrow(data_test),ncol=iter+3))
  colnames(result_for_test_detail)<-c(paste('iter',1:iter,sep='_'),'city','vote','check')
  row.names(result_for_test_detail)<-row.names(data_test)
  for(i in 1:nrow(result_for_test_detail)){
    result_for_test_detail[i,'city']<-tiqu(row.names(result_for_test_detail)[i])
  }
  for(i in 1:iter){
    ## Cross-Validation
    set.seed(i)
    k<-sampling(data_training,k1)
    data_train<-data_training[-k,]
    data_validation<-data_training[k,]
    lda_fit<-lda(factor(city)~.,data=data_train)
    lda_pred<-predict(lda_fit,data_validation)
    lda_predict<-as.data.frame(cbind(row.names(data_training)[k],as.vector(lda_pred$class)))
    colnames(lda_predict)<-c('sample','predict')
    lda_predict$sample<-as.character(lda_predict$sample)
    lda_predict$predict<-as.character(lda_predict$predict)
    lda_predict$city<-0
    for(j in 1:nrow(lda_predict)){
      lda_predict[j,3]<-tiqu(lda_predict[j,1])
    }
    rm(j)
    lda_predict$check<-0
    for(j in 1:nrow(lda_predict)){
      if(lda_predict[j,'predict']==lda_predict[j,'city']){lda_predict[j,'check']<-1}
    }
    rm(j)
    lda_predict<-lda_predict[,c(1,3,2,4)]
    lda_predict_city<-as.data.frame(matrix(0,nrow=length(city_all),ncol=2))
    colnames(lda_predict_city)<-c('city','error_rate')
    for(j in 1:length(city_all)){
      lda_predict_city[j,1]<-city_all[j]
      lda_predict_city[j,2]<-1-mean(lda_predict[which(lda_predict$city==city_all[j]),'check'])
    }
    rm(j)
    lda_predict_city<-lda_predict_city[order(lda_predict_city$error_rate,decreasing = F),]
    error_rate<-1-mean(lda_predict$check)
    result_for_1_iter_CV<-list(error_rate=error_rate,
                               lda_predict=lda_predict,
                               lda_predict_city=lda_predict_city)
    result_for_CV[[i]]<-result_for_1_iter_CV
    
    ## Test prediction
    lda_pred_test<-predict(lda_fit,data_test)
    result_for_test_detail[,i]<-as.vector(lda_pred_test)$class
  }
  ## vote result
  result_for_test_detail$vote<-apply(result_for_test_detail[,1:iter],1,getmode)
  for(i in 1:nrow(result_for_test_detail)){
    if(result_for_test_detail[i,'vote']==result_for_test_detail[i,'city']){result_for_test_detail[i,'check']<-1}
  }
  result_for_test<-result_for_test_detail[,c('city','vote','check')]
  time_end<-Sys.time()
  time_diff<-time_end-time_start
  print(time_diff)
  result_output<-list(result_for_CV=result_for_CV,
                      result_for_test_detail=result_for_test_detail,
                      result_for_test=result_for_test,
                      time=time_diff)
  return(result_output)
}

#### function: ml ####
ml_f<-function(data_input,seed,k_validation,k_test,iter,ml){
  data_normalized<-normalization_f(data_input)
  set.seed(seed)
  test_ind<-sampling(data_normalized,k_test)
  data_test<-data_normalized[test_ind,]
  data_training<-data_normalized[-test_ind,]
  output_result<-list()
  if('rf'%in%ml){
    ## random forest
    rf_result<-rf_k_fold(data_training,data_test,k_validation,iter)
    output_result[[grep('rf',ml)]]<-rf_result
  }
  if('svm'%in%ml){
    ## svm
    svm_result<-svm_k_fold(data_training,data_test,k_validation,iter)
    output_result[[grep('svm',ml)]]<-svm_result
  }
  if('lda'%in%ml){
    ## lda
    lda_result<-lda_k_fold(data_training,data_test,k_validation,iter)
    output_result[[grep('lda',ml)]]<-lda_result
  }
  if('deep'%in%ml){
    ## deep learning
    deep_result<-deep_learning_k_fold(data_training,data_test,k_validation,iter)
    output_result[[grep('deep',ml)]]<-deep_result
  }
  return(output_result)
}

#### Function: error rate ####
error_rate<-function(data,ml){
  ml_type<-c('rf','svm','deep')
  error_rate_table<-as.data.frame(matrix(NA,ncol=length(ml),nrow=1))
  error_rate_city_table<-as.data.frame(matrix(NA,ncol=length(ml)+1,nrow=28))
  colnames(error_rate_table)<-ml
  colnames(error_rate_city_table)<-c('city',ml)
  error_rate_city_table$city<-unique(data[[1]]$result_for_test$city)
  for(i in 1:length(ml)){
    index<-which(ml_type==ml[i])
    error_rate_table[1,i]<-1-mean(data[[index]]$result_for_test$check)
    bb<-data[[index]]$result_for_test %>% group_by(city) %>% summarise(error_rate=1-mean(check))
    error_rate_city_table[,i+1]<-bb$error_rate
  }
  result_output<-list(error_rate_table=error_rate_table,
                      error_rate_city_table=error_rate_city_table)
}


#### Function: validation error rate ####
validation_error_rate<-function(data_list,iter){
  validation_error_rate_table<-as.data.frame(matrix(NA,nrow=iter+1,ncol=3))
  colnames(validation_error_rate_table)<-c('rf','svm','deep')
  for(i in 1:iter){
    validation_error_rate_table[i,1]<-data_list[[1]]$result_for_CV[[iter]]$error_rate
    validation_error_rate_table[i,2]<-data_list[[2]]$result_for_CV[[iter]]$error_rate
    validation_error_rate_table[i,3]<-1-data_list[[3]]$all_history[[iter]]$metrics$val_accuracy[200]
  }
  validation_error_rate_table[iter+1,]<-apply(validation_error_rate_table[1:iter,],2,mean)
  return(validation_error_rate_table)
}
#### validation error rate result ####
validation_error_rate(result_10_123,100)[101,]
validation_error_rate(result_30_123,100)[101,]
validation_error_rate(result_50_123,100)[101,]
validation_error_rate(result_70_123,100)[101,]
validation_error_rate(result_100_123,100)[101,]

#### city error rate plot ####
data_30_for_plot<-error_rate(result_30_123,c('rf','svm','deep'))$error_rate_city_table
data_30_for_plot_long<-gather(data_30_for_plot,method,error_rate,c('rf','svm','deep'))
temp_data<-data_30_for_plot_long %>%
  group_by(city) %>%
  summarise(mean=mean(error_rate))
data_30_for_plot_long<-inner_join(data_30_for_plot_long,temp_data,by='city')

data_30_for_plot_long$method<-factor(data_30_for_plot_long$method, levels=c('rf','svm','deep'), labels=c('RF','SVM','MLP')) 
plot_error_city<-ggplot(data=data_30_for_plot_long,aes(x=reorder(city,mean),y=error_rate,fill=city))+
  facet_grid(.~method)+
  geom_bar(stat='identity')+
  xlab('city')+
  ylab('error rate')+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))

print(plot_error_city)
  
#### data save for PCoA ####
save(list=c('data_used_normalized_12',
            'all_city_samples_table',
            'data_used_12',
            'normalization_f'
),
file='Main dataset 7 data for PCoA.RData')

#### data save ####
save(list=c('all_city_samples_table',
            'data_test_12',
            'data_training_12',
            'data_used_12',
            'data_used_normalized_12',
            'all_city_samples_filtered',
            'city_all',
            'test_ind',
            'deep_learning_k_fold',
            'getmode',
            'ml_f',
            'normalization_f',
            'rf_k_fold',
            'sampling',
            'svm_k_fold',
            'tiqu'
            ),file='Main dataset 6 data for multi-output model.RData')

save('result_deep_12',file='result_deep_12.RData')

## 12 ##
data_used_12<-features_0.6_otu_table[features_used_12,]
data_used_normalized_12<-normalization_f(data_used_12)
data_test_12<-data_used_normalized_12[test_ind,]
data_training_12<-data_used_normalized_12[-test_ind,]
result_rf_12<-rf_k_fold(data_training_12,data_test_12,5,100)
result_svm_12<-svm_k_fold(data_training_12,data_test_12,5,100)
result_deep_12_test<-deep_learning_k_fold(data_training_12,data_test_12,5,100)
# result_ml_12<-ml_f(data_used,123,5,5,100,c('rf','svm','deep'))

#### vote1 and vote2 ####
data_deep_12<-result_deep_12$result_for_test_detail
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
data_deep_12$vote_number<-apply(data_deep_12[,1:100],1,get_mode_number)
data_deep_12$vote_2<-apply(data_deep_12[,1:100],1,get_second_vote)
data_deep_12$vote_2_number<-apply(data_deep_12[,1:100],1,get_second_vote_number)
data_deep_12$vote_2[20]<-'ARN'
data_deep_12$vote_2[41]<-'DEN_17  '
data_deep_12$vote_2[109]<-'IEV'

data_deep_12_result<-data_deep_12[,101:106]
data_deep_12_result$check_2<-0
data_deep_12_result$check_2[which(data_deep_12_result$city==data_deep_12_result$vote_2)]<-1
data_deep_12_result$check_1_2<-data_deep_12_result$check+data_deep_12_result$check_2

temp_data<-data_deep_12_result %>%
  group_by(city) %>%
  summarise(mean_1=1-mean(check))
data_deep_12_result<-inner_join(data_deep_12_result,temp_data,by='city')

temp_data<-data_deep_12_result %>%
  group_by(city) %>%
  summarise(mean_1_2=1-mean(check_1_2))
data_deep_12_result<-inner_join(data_deep_12_result,temp_data,by='city')

#### Figure 1 ####
plot_error_city_1<-ggplot(data=unique(data_deep_12_result[c('city','mean_1')]),aes(x=reorder(city,mean_1),y=mean_1,fill=city))+
  geom_bar(stat='identity')+
  xlab('city')+
  ylab('error rate')+
  ggtitle('Main dataset vote 1')+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))

plot_error_city_1_2<-ggplot(data=unique(data_deep_12_result[c('city','mean_1_2')]),aes(x=reorder(city,mean_1_2),y=mean_1_2,fill=city))+
  geom_bar(stat='identity')+
  xlab('city')+
  ylab('error rate')+
  ggtitle('Main dataset vote 1_2')+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))

print(plot_error_city_1)
print(plot_error_city_1_2)
