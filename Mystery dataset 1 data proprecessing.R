#### read data ####
setwd('E:/University of Florida/My research/CAMDA 2020/mystery/data')
mystery_samples_list<-read.table('files2.txt')
mystery_samples_list$mystery_ID<-paste0('mystery_',sprintf("%03d", 1:nrow(mystery_samples_list)))
mystery_samples<-sprintf("%03d", 1:nrow(mystery_samples_list))

for(i in 1:nrow(mystery_samples_list)){
  eval(parse(text=paste0(mystery_samples_list[i,2],'<-read.delim(\'',mystery_samples_list[i,1],'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',mystery_samples[i],'<-read.delim(\'',mystery_samples_list[i,1],'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

#### for large sample ####
for(j in 1:nrow(mystery_samples_list)){
  eval(parse(text=paste0('mystery','_',mystery_samples[j],'<-f_large_sample(','\'mystery\',\'',mystery_samples[j],'\')')))
}

for(j in 1:nrow(mystery_samples_list)){
  otu_info<-eval(parse(text=paste0('otu_',mystery_samples[j])))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',mystery_samples[j],'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

#### data proprecessing ####
for(i in 1:nrow(mystery_samples_list)){
  eval(parse(text=paste0('otu_',mystery_samples[i],'[,2]<-gsub(\'[][]\',\'\',otu_',mystery_samples[i],'[,2])')))
}

for(i in 1:nrow(mystery_samples_list)){
  eval(parse(text=paste0('mystery','_',mystery_samples[i],'<-OTU_generate(','\'mystery\',\'',mystery_samples[i],'\')')))
}

for(i in 1:nrow(mystery_samples_list)){
  eval(parse(text=paste0('rm(otu_',mystery_samples[i],')')))
}

mystery_species<-NULL
for(i in mystery_samples){
  bb<-eval(parse(text=paste0('unique(mystery_',i,'$Species)')))
  mystery_species<-union(mystery_species,bb)
}
mystery_species<-na.omit(mystery_species)

mystery_family<-NULL
for(i in mystery_samples){
  bb<-eval(parse(text=paste0('unique(mystery_',i,'$Family)')))
  mystery_family<-union(mystery_family,bb)
}
mystery_family<-na.omit(mystery_family)

mystery_order<-NULL
for(i in mystery_samples){
  bb<-eval(parse(text=paste0('unique(mystery_',i,'$Order)')))
  mystery_order<-union(mystery_order,bb)
}
mystery_order<-na.omit(mystery_order)

rm(i,j,bb)
#### otu table ####
mystery_otu_species<-otu_table_generate('mystery','species',mystery_samples)
mystery_otu_family<-otu_table_generate('mystery','family',mystery_samples)
mystery_otu_order<-otu_table_generate('mystery','order',mystery_samples)

setwd('E:/University of Florida/My research/CAMDA 2020/data')
save(list=c('mystery_otu_order',
            'mystery_otu_family',
            'mystery_otu_species',
            'mystery_samples_list',
            'mystery_samples'
),
file='Mystery data for analysis')


#### name of the samples ####
setwd('E:/University of Florida/My research/CAMDA 2020/mystery')
save(list=c('mystery_samples_list'),
     file='mystery_samples_list.RData')



