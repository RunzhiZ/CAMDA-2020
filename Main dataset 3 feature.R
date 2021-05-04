load('Main dataset 2 data ready for analysis.RData')

#### Function: remove the samples with low quality ####
rm_median_p<-function(data_input,p){
  data_used<-apply(data_input,2,sum)
  data_median<-median(data_used)
  samples_filtered<-names(data_used)[which(data_used>=data_median*p)]
  return(samples_filtered)
}

all_city_samples_filtered<-NULL
all_city_samples_filtered_table<-as.data.frame(matrix(NA,nrow=28,ncol=3))
colnames(all_city_samples_filtered_table)<-c('city','samples','samples_filtered')
for(i in 1:length(city_all)){
  aa<-rm_median_p(get(paste0(city_all[i],'_otu_order')),0.2)
  all_city_samples_filtered_table[i,1]<-city_all[i]
  bb<-get(paste0(city_all[i],'_otu_order'))
  all_city_samples_filtered_table[i,2]<-ncol(bb)
  all_city_samples_filtered_table[i,3]<-length(aa)
  all_city_samples_filtered<-c(all_city_samples_filtered,aa)
  for(j in taxa)
  eval(parse(text=paste0(city_all[i],'_otu_',j,'_filtered','<-',city_all[i],'_otu_',j,'[,aa]')))
  rm(aa,bb)
}

#### calculate the ubiquity for each city ####
list_ubiquity_order<-list()
list_ubiquity_family<-list()
#list_ubiquity_genus<-list()
list_ubiquity_species<-list()
taxa<-c('species','family','order')
for(i in city_all){
  for(j in taxa){
    aa<-apply((get(paste(i,'otu',j,'filtered',sep='_'))!=0),1,sum)/ncol(get(paste(i,'otu',j,'filtered',sep='_')))
    eval(parse(text=paste0('list_ubiquity_',j,'$',i,'<-aa')))
    rm(aa)
  }
}

#### ubiquity table for different taxa ####
#### species
all_species_ubiquity_table<-as.data.frame(matrix(0,nrow=length(all_species),ncol=length(city_all)))
row.names(all_species_ubiquity_table)<-all_species
colnames(all_species_ubiquity_table)<-city_all
#### genus
#all_genus_ubiquity_table<-as.data.frame(matrix(0,nrow=length(all_genus),ncol=length(city_all)))
#row.names(all_genus_ubiquity_table)<-all_genus
#colnames(all_genus_ubiquity_table)<-city_all
#### family
all_family_ubiquity_table<-as.data.frame(matrix(0,nrow=length(all_family),ncol=length(city_all)))
row.names(all_family_ubiquity_table)<-all_family
colnames(all_family_ubiquity_table)<-city_all
#### order
all_order_ubiquity_table<-as.data.frame(matrix(0,nrow=length(all_order),ncol=length(city_all)))
row.names(all_order_ubiquity_table)<-all_order
colnames(all_order_ubiquity_table)<-city_all

for(i in city_all){
  for(j in taxa){
    eval(parse(text=paste0('all_',j,'_ubiquity_table[names(list_ubiquity_',j,'$',i,'),\'',i,'\']<-list_ubiquity_',j,'$',i)))
  }
}

#### Function: Obtain the features based on their ubiquity ####
## Calculate the number of cities with ubiquity > 0.xx
length_0.9<-function(data){
  return(length(which(data>0.9)))
}
length_0.8<-function(data){
  return(length(which(data>0.8)))
}
length_0.7<-function(data){
  return(length(which(data>0.7)))
}
length_0.6<-function(data){
  return(length(which(data>0.6)))
}
length_0.5<-function(data){
  return(length(which(data>0.5)))
}

#### calculate the number of  cities using different ubiquity cutoff ####
all_species_ubiquity_table$ubiquity_0.5<-apply(all_species_ubiquity_table[,1:28],1,length_0.5)
all_species_ubiquity_table$ubiquity_0.6<-apply(all_species_ubiquity_table[,1:28],1,length_0.6)
all_species_ubiquity_table$ubiquity_0.7<-apply(all_species_ubiquity_table[,1:28],1,length_0.7)
all_species_ubiquity_table$ubiquity_0.8<-apply(all_species_ubiquity_table[,1:28],1,length_0.8)
all_species_ubiquity_table$ubiquity_0.9<-apply(all_species_ubiquity_table[,1:28],1,length_0.9)

#all_genus_ubiquity_table$ubiquity_0.5<-apply(all_genus_ubiquity_table[,1:28],1,length_0.5)
#all_genus_ubiquity_table$ubiquity_0.6<-apply(all_genus_ubiquity_table[,1:28],1,length_0.6)
#all_genus_ubiquity_table$ubiquity_0.7<-apply(all_genus_ubiquity_table[,1:28],1,length_0.7)
#all_genus_ubiquity_table$ubiquity_0.8<-apply(all_genus_ubiquity_table[,1:28],1,length_0.8)
#all_genus_ubiquity_table$ubiquity_0.9<-apply(all_genus_ubiquity_table[,1:28],1,length_0.9)

all_family_ubiquity_table$ubiquity_0.5<-apply(all_family_ubiquity_table[,1:28],1,length_0.5)
all_family_ubiquity_table$ubiquity_0.6<-apply(all_family_ubiquity_table[,1:28],1,length_0.6)
all_family_ubiquity_table$ubiquity_0.7<-apply(all_family_ubiquity_table[,1:28],1,length_0.7)
all_family_ubiquity_table$ubiquity_0.8<-apply(all_family_ubiquity_table[,1:28],1,length_0.8)
all_family_ubiquity_table$ubiquity_0.9<-apply(all_family_ubiquity_table[,1:28],1,length_0.9)

all_order_ubiquity_table$ubiquity_0.5<-apply(all_order_ubiquity_table[,1:28],1,length_0.5)
all_order_ubiquity_table$ubiquity_0.6<-apply(all_order_ubiquity_table[,1:28],1,length_0.6)
all_order_ubiquity_table$ubiquity_0.7<-apply(all_order_ubiquity_table[,1:28],1,length_0.7)
all_order_ubiquity_table$ubiquity_0.8<-apply(all_order_ubiquity_table[,1:28],1,length_0.8)
all_order_ubiquity_table$ubiquity_0.9<-apply(all_order_ubiquity_table[,1:28],1,length_0.9)

#### all features table ####
all_species_table<-as.data.frame(matrix(0,nrow=length(all_species),ncol=length(all_city_samples_filtered)))
row.names(all_species_table)<-all_species
colnames(all_species_table)<-all_city_samples_filtered
#all_genus_table<-as.data.frame(matrix(0,nrow=length(all_genus),ncol=length(all_city_samples_filtered)))
#row.names(all_genus_table)<-all_genus
#colnames(all_genus_table)<-all_city_samples_filtered
all_family_table<-as.data.frame(matrix(0,nrow=length(all_family),ncol=length(all_city_samples_filtered)))
row.names(all_family_table)<-all_family
colnames(all_family_table)<-all_city_samples_filtered
all_order_table<-as.data.frame(matrix(0,nrow=length(all_order),ncol=length(all_city_samples_filtered)))
row.names(all_order_table)<-all_order
colnames(all_order_table)<-all_city_samples_filtered

for(i in city_all){
  for(j in taxa){
    eval(parse(text=paste0('all_',j,'_table[row.names(',i,'_otu_',j,'_filtered),colnames(',i,'_otu_',j,'_filtered)]<-',i,'_otu_',j,'_filtered')))
  }
}

#### Function: get the group of the samples####
tiqu<-function(str){
  substr(str,1,min(which(strsplit(str,split='')[[1]]=='0'))-2)
}
#### all city samples table ####
all_city_samples_table<-as.data.frame(matrix(NA,nrow=length(all_city_samples_filtered),ncol=3))
all_city_samples_table[,1]<-all_city_samples_filtered
for(i in 1:nrow(all_city_samples_table)){
  all_city_samples_table[i,2]<-tiqu(all_city_samples_table[i,1])
}
all_city_samples_table[,3]<-substr(all_city_samples_table[,2],1,3)
colnames(all_city_samples_table)<-c('sample','group_1','group_2')

#### Number of species/family/order in the filtered data ####
length(which(apply(all_order_table,1,sum)!=0))
length(which(apply(all_family_table,1,sum)!=0))
length(which(apply(all_species_table,1,sum)!=0))

#### save data ####
save.image('Main dataset 3 feature.RData')

#### data for machine learning ####
save(list=c('all_city_samples_table',
            'all_family_table',
            'all_family_ubiquity_table',
            'all_order_table',
            'all_order_ubiquity_table',
            'all_species_table',
            'all_species_ubiquity_table',
            'all_city_samples_filtered',
            'city_all',
            'taxa',
            'tiqu'
            ),
     file='Main dataset 4 data for feature selection and machine learning.RData')


#### data overview ####
data_overview<-data.frame(matrix(0,nrow=length(city_all),ncol=6))
colnames(data_overview)<-c('city','number_of_samples','species','genus','family','order')
data_overview$city<-city_all
for(i in 1:length(city_all)){
  data_overview$number_of_samples[i]<-length(grep(city_all[i],all_city_samples))
  data_overview$species[i]<-nrow(get(paste0(city_all[i],'_otu_','species')))
  data_overview$genus[i]<-nrow(get(paste0(city_all[i],'_otu_','genus')))
  data_overview$family[i]<-nrow(get(paste0(city_all[i],'_otu_','family')))
  data_overview$order[i]<-nrow(get(paste0(city_all[i],'_otu_','order')))
}
write.csv(data_overview,'data_overview.csv')




