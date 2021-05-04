load('Main dataset 1 Generate OTU tables for all samples.RData')
#### Function: generate otu table ####
otu_table_generate<-function(city,taxa,samples){
  n<-eval(parse(text=paste0('length(',city,'_',taxa,')')))
  otu_table<-as.data.frame(matrix(nrow=n,ncol=length(samples)))
  city_taxa<-eval(parse(text=paste0(city,'_',taxa)))
  row.names(otu_table)<-city_taxa
  for(i in 1:length(samples)){
    colnames(otu_table)[i]<-paste0(city,'_',samples[i])
    single_city<-eval(parse(text=paste0(city,'_',samples[i])))
    if(taxa=='species'){
      for(j in city_taxa){
        if(length(which(single_city$Species %in% j))==0){otu_table[j,i]<-0}
        else{otu_table[j,i]<-sum(single_city[which(single_city$Species==j),'count'])}
      }
    }
    if(taxa=='genus'){
      for(j in city_taxa){
        if(length(which(single_city$Genus %in% j))==0){otu_table[j,i]<-0}
        else{otu_table[j,i]<-sum(single_city[which(single_city$Genus==j),'count'])}
      }
    }
    if(taxa=='family'){
      for(j in city_taxa){
        if(length(which(single_city$Family %in% j))==0){otu_table[j,i]<-0}
        else{otu_table[j,i]<-sum(single_city[which(single_city$Family==j),'count'])}
      }
    }
    if(taxa=='order'){
      for(j in city_taxa){
        if(length(which(single_city$Order %in% j))==0){otu_table[j,i]<-0}
        else{otu_table[j,i]<-sum(single_city[which(single_city$Order==j),'count'])}
      }
    }
  }
  return(otu_table)
}

#### Obtain the common species, families and orders ####
l_species<-list(as.vector(ARN_species),
                as.vector(BCN_species),
                as.vector(BER_species),
                as.vector(DEN_16_species),
                as.vector(DEN_17_species),
                as.vector(DOH_16_species),
                as.vector(DOH_17_species),
                as.vector(FAI_species),
                as.vector(HKG_species),
                as.vector(ICN_species),
                as.vector(IEV_species),
                as.vector(ILR_16_species),
                as.vector(ILR_17_species),
                as.vector(KUL_species),
                as.vector(LCY_species),
                as.vector(LIS_species),
                as.vector(NYC_16_species),
                as.vector(NYC_17_species),
                as.vector(OFF_species),
                as.vector(SAO_species),
                as.vector(SCL_species),
                as.vector(SDJ_species),
                as.vector(SFO_species),
                as.vector(SGP_species),
                as.vector(TPE_species),
                as.vector(TYO_16_species),
                as.vector(TYO_17_species),
                as.vector(ZRH_species))

l_genus<-list(as.vector(ARN_genus),
                as.vector(BCN_genus),
                as.vector(BER_genus),
                as.vector(DEN_16_genus),
                as.vector(DEN_17_genus),
                as.vector(DOH_16_genus),
                as.vector(DOH_17_genus),
                as.vector(FAI_genus),
                as.vector(HKG_genus),
                as.vector(ICN_genus),
                as.vector(IEV_genus),
                as.vector(ILR_16_genus),
                as.vector(ILR_17_genus),
                as.vector(KUL_genus),
                as.vector(LCY_genus),
                as.vector(LIS_genus),
                as.vector(NYC_16_genus),
                as.vector(NYC_17_genus),
                as.vector(OFF_genus),
                as.vector(SAO_genus),
                as.vector(SCL_genus),
                as.vector(SDJ_genus),
                as.vector(SFO_genus),
                as.vector(SGP_genus),
                as.vector(TPE_genus),
                as.vector(TYO_16_genus),
                as.vector(TYO_17_genus),
                as.vector(ZRH_genus))

l_family<-list(as.vector(ARN_family),
                as.vector(BCN_family),
                as.vector(BER_family),
                as.vector(DEN_16_family),
                as.vector(DEN_17_family),
                as.vector(DOH_16_family),
                as.vector(DOH_17_family),
                as.vector(FAI_family),
                as.vector(HKG_family),
                as.vector(ICN_family),
                as.vector(IEV_family),
                as.vector(ILR_16_family),
                as.vector(ILR_17_family),
                as.vector(KUL_family),
                as.vector(LCY_family),
                as.vector(LIS_family),
                as.vector(NYC_16_family),
                as.vector(NYC_17_family),
                as.vector(OFF_family),
                as.vector(SAO_family),
                as.vector(SCL_family),
                as.vector(SDJ_family),
                as.vector(SFO_family),
                as.vector(SGP_family),
                as.vector(TPE_family),
                as.vector(TYO_16_family),
                as.vector(TYO_17_family),
                as.vector(ZRH_family))

l_order<-list(as.vector(ARN_order),
                as.vector(BCN_order),
                as.vector(BER_order),
                as.vector(DEN_16_order),
                as.vector(DEN_17_order),
                as.vector(DOH_16_order),
                as.vector(DOH_17_order),
                as.vector(FAI_order),
                as.vector(HKG_order),
                as.vector(ICN_order),
                as.vector(IEV_order),
                as.vector(ILR_16_order),
                as.vector(ILR_17_order),
                as.vector(KUL_order),
                as.vector(LCY_order),
                as.vector(LIS_order),
                as.vector(NYC_16_order),
                as.vector(NYC_17_order),
                as.vector(OFF_order),
                as.vector(SAO_order),
                as.vector(SCL_order),
                as.vector(SDJ_order),
                as.vector(SFO_order),
                as.vector(SGP_order),
                as.vector(TPE_order),
                as.vector(TYO_16_order),
                as.vector(TYO_17_order),
                as.vector(ZRH_order))

common_species<-Reduce(intersect,l_species) ## 7
common_genus<-Reduce(intersect,l_genus) ## 6
common_family<-Reduce(intersect,l_family) ## 9
common_order<-Reduce(intersect,l_order) ## 12
all_species<-Reduce(union,l_species) ## 1394
all_genus<-Reduce(union,l_genus) ## 863
all_family<-Reduce(union,l_family) ## 288
all_order<-Reduce(union,l_order) ## 193

#### calculate the number of feature ####
city_all<-c('ARN','BCN','BER','DEN_16','DEN_17','DOH_16','DOH_17','FAI','HKG','ICN','IEV','ILR_16','ILR_17','KUL','LCY',
            'LIS','NYC_16','NYC_17','OFF','SAO','SCL','SDJ','SFO','SGP','TPE','TYO_16','TYO_17','ZRH')
city_species<-data.frame()
for(i in city_all){
  city_species[which(city_all==i),1]<-i
  city_species[which(city_all==i),2]<-eval(parse(text=paste0('length(',i,'_species)')))
}
colnames(city_species)<-c('city','number_species')

city_genus<-data.frame()
for(i in city_all){
  city_genus[which(city_all==i),1]<-i
  city_genus[which(city_all==i),2]<-eval(parse(text=paste0('length(',i,'_genus)')))
}
colnames(city_genus)<-c('city','number_genus')

city_family<-data.frame()
for(i in city_all){
  city_family[which(city_all==i),1]<-i
  city_family[which(city_all==i),2]<-eval(parse(text=paste0('length(',i,'_family)')))
}
colnames(city_family)<-c('city','number_family')

city_order<-data.frame()
for(i in city_all){
  city_order[which(city_all==i),1]<-i
  city_order[which(city_all==i),2]<-eval(parse(text=paste0('length(',i,'_order)')))
}
colnames(city_order)<-c('city','number_order')

city_feature<-cbind(city_species,city_genus$number_genus,city_family$number_family,city_order$number_order)
colnames(city_feature)[c(3:5)]<-c('number_genus','number_family','number_order')

#### Species ####
##
##
##
##
for(i in city_all){
  temp_table<-otu_table_generate(i,'species',get(paste0(i,'_','samples')))
  eval(parse(text=paste0(i,'_','otu','_','species','<-','temp_table')))
  rm(temp_table)
}

#### Genus ####
##
##
##
##
for(i in city_all){
  temp_table<-otu_table_generate(i,'genus',get(paste0(i,'_','samples')))
  eval(parse(text=paste0(i,'_','otu','_','genus','<-','temp_table')))
  rm(temp_table)
}

#### Family ####
##
##
##
##
for(i in city_all){
  temp_table<-otu_table_generate(i,'family',get(paste0(i,'_','samples')))
  eval(parse(text=paste0(i,'_','otu','_','family','<-','temp_table')))
  rm(temp_table)
}

#### Order ####
##
##
##
##
for(i in city_all){
  temp_table<-otu_table_generate(i,'order',get(paste0(i,'_','samples')))
  eval(parse(text=paste0(i,'_','otu','_','order','<-','temp_table')))
  rm(temp_table)
}

#### all_city_samples ####
all_city_samples<-NULL
for(i in city_all){
  aa<-get(paste0(i,'_','samples'))
  bb<-paste0(i,'_',aa)
  all_city_samples<-c(all_city_samples,bb)
}

## otu all species ##
otu_all_species<-as.data.frame(matrix(0,nrow=length(all_species),ncol=length(all_city_samples)))
row.names(otu_all_species)<-all_species
colnames(otu_all_species)<-all_city_samples
for(i in city_all){
  otu_all_species[row.names(get(paste0(i,'_otu_','species'))),colnames(get(paste0(i,'_otu_','species')))]<-get(paste0(i,'_otu_','species'))
}

## otu all genus ##
otu_all_genus<-as.data.frame(matrix(0,nrow=length(all_genus),ncol=length(all_city_samples)))
row.names(otu_all_genus)<-all_genus
colnames(otu_all_genus)<-all_city_samples
for(i in city_all){
  otu_all_genus[row.names(get(paste0(i,'_otu_','genus'))),colnames(get(paste0(i,'_otu_','genus')))]<-get(paste0(i,'_otu_','genus'))
}

## otu all family ##
otu_all_family<-as.data.frame(matrix(0,nrow=length(all_family),ncol=length(all_city_samples)))
row.names(otu_all_family)<-all_family
colnames(otu_all_family)<-all_city_samples
for(i in city_all){
  otu_all_family[row.names(get(paste0(i,'_otu_','family'))),colnames(get(paste0(i,'_otu_','family')))]<-get(paste0(i,'_otu_','family'))
}

## otu all order ##
otu_all_order<-as.data.frame(matrix(0,nrow=length(all_order),ncol=length(all_city_samples)))
row.names(otu_all_order)<-all_order
colnames(otu_all_order)<-all_city_samples
for(i in city_all){
  otu_all_order[row.names(get(paste0(i,'_otu_','order'))),colnames(get(paste0(i,'_otu_','order')))]<-get(paste0(i,'_otu_','order'))
}

#### remove unnecessary variables ####
#### remove unnecessary variables ####
#### remove unnecessary variables ####
for(i in city_all){
  for(j in get(paste0(i,'_samples'))){
    eval(parse(text=paste0('rm(',i,'_',j,')')))
  }
}
taxa<-c('species','genus','family','order')
for(i in city_all){
  for(j in taxa){
    eval(parse(text=paste0('rm(',i,'_',j,')')))
  }
}
for(i in city_all){
  eval(parse(text=paste0('rm(',i,'_','samples',')')))
}
for(i in city_all){
  eval(parse(text=paste0('rm(',i,'_otu_','genus',')')))
}

#### save data ####
save.image('Main dataset 2 data ready for analysis.RData')