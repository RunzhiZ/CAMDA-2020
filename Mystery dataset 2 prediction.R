library(limma)
library(edgeR)
library(dplyr)
library(keras)

feature_species<-intersect(row.names(data_used_12),row.names(mystery_otu_species))
feature_family<-intersect(row.names(data_used_12),row.names(mystery_otu_family))
feature_order<-intersect(row.names(data_used_12),row.names(mystery_otu_order))
mystery_otu_feature<-rbind(mystery_otu_species[feature_species,],
                           mystery_otu_family[feature_family,],
                           mystery_otu_order[feature_order,])
mystery_otu_feature_filtered<-mystery_otu_feature[,-which(colSums(mystery_otu_feature)<=median(colSums(mystery_otu_feature))/5)]
mix_otu_feature<-cbind(data_used_12,mystery_otu_feature_filtered)

#### normalization ####
d0<-DGEList(mix_otu_feature)
group<-c(all_city_samples_table[which(all_city_samples_table[,1] %in% colnames(mix_otu_feature)),2],rep('mystery',ncol(mystery_otu_feature_filtered)))
mm<-model.matrix(~0 + group)
y<-voom(d0, mm, plot = F)$E
otu_normalized<-as.data.frame(t(y))
otu_normalized$city<-group
names(otu_normalized)<-make.names(names(otu_normalized))
rm(d0,y,mm,group)

#### deep learning ####
id_table<-data.frame(city=unique(otu_normalized$city),id=0:28)
otu_normalized_id<-inner_join(otu_normalized,id_table,by=c('city'='city'))
row.names(otu_normalized_id)<-row.names(otu_normalized)

result_table_all<-data.frame()
for(i in 1:100){
  set.seed(i)
  mix_test_ind<-sample(993:1103,22,replace = F)
  mix_training<-otu_normalized_id[-mix_test_ind,]
  mix_test<-otu_normalized_id[mix_test_ind,]
  mix_training_Labels <- to_categorical(mix_training[,'id'])
  #### deep learning
  model <- keras_model_sequential() 
  
  # Add layers to the model
  model %>% 
    layer_dense(units = 128, activation = 'relu', input_shape = c(44)) %>% 
    layer_dense(units = 128, activation = 'relu') %>%
    layer_dropout(rate=0.5) %>%
    layer_dense(units = 128, activation = 'relu') %>%
    layer_dropout(rate=0.5) %>%
    layer_dense(units = 29, activation = 'softmax')
  
  # Compile the model
  model %>% compile(
    loss = 'categorical_crossentropy',
    optimizer = 'adam',
    metrics = 'accuracy'
  )
  
  # Store the fitting history in `history` 
  history <- model %>% fit(
    as.matrix(mix_training[,1:44]), 
    mix_training_Labels, 
    epochs = 500,
    batch_size = 32, 
    verbose=0
  )
  # Predict the classes for the test data
  classes <- model %>% predict_classes(as.matrix(mix_test[,1:44]), batch_size = 32)
  result_table<-cbind(row.names(mix_test),classes)
  result_table_all<-rbind(result_table_all,result_table)
}

#### read true label ####
mystery_true<-read.csv('CAMDA_mysterySample_City.csv',header=T)
mystery_samples_list<-inner_join(mystery_samples_list,mystery_true,by=c('V1'='New_name'))
mystery_samples_list$city_code<-mystery_samples_list$City
mystery_samples_list[which(mystery_samples_list$City=='hong_kong'),'city_code']<-'HKG'
mystery_samples_list[which(mystery_samples_list$City=='tokyo'),'city_code']<-'TYO'
mystery_samples_list[which(mystery_samples_list$City=='taipei'),'city_code']<-'TPE'
mystery_samples_list[which(mystery_samples_list$City=='zurich'),'city_code']<-'ZRH'
mystery_samples_list[which(mystery_samples_list$City=='kyiv'),'city_code']<-'IEV'

colnames(result_table_all)[1]<-'sample'
mystery_table<-as.data.frame(matrix(NA,nrow=111,ncol=5))
colnames(mystery_table)<-c('sample','vote_1','vote_1_number','vote_2','vote_2_number')
mystery_table$sample<-colnames(mystery_otu_feature_filtered)
for(i in 1:nrow(mystery_table)){
  sample_table<-table(result_table_all[which(result_table_all$sample==mystery_table$sample[i]),][,2])
  sample_table_order<-sample_table[order(sample_table,decreasing = T)]
  mystery_table$vote_1[i]<-names(sample_table_order)[1]
  mystery_table$vote_2[i]<-names(sample_table_order)[2]
  mystery_table$vote_1_number[i]<-sample_table_order[1]/sum(sample_table)
  mystery_table$vote_2_number[i]<-sample_table_order[2]/sum(sample_table)
}
mystery_table$vote_1<-as.numeric(mystery_table$vote_1)
mystery_table$vote_2<-as.numeric(mystery_table$vote_2)
mystery_table<-inner_join(mystery_table,id_table,by=c('vote_1'='id'))
mystery_table<-inner_join(mystery_table,id_table,by=c('vote_2'='id'))
mystery_table<-inner_join(mystery_table,mystery_samples_list[,c(2,4)],by=c('sample'='mystery_ID'))
colnames(mystery_table)[6:8]<-c('pred_1','pred_2','true')


