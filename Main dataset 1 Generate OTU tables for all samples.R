#### Function ####
## For large samples
f_large_sample<-function(city,city_sample){
  if(ncol(get(paste0(city,'_',city_sample)))==3){
    bb<-data.frame()
    otu_table<-eval(parse(text=paste0(city,'_',city_sample,'[c(-1,-2),]')))
    uniq<-unique(otu_table$V2)
    uniq<-uniq[uniq!='']
    for(i in uniq){
      bb[which(uniq==i),1]<-i
      bb[which(uniq==i),2]<-sum(otu_table[which(otu_table$V2==i),'V3'])
    }
    return(bb)
  }else(return(get(paste0(city,'_',city_sample))))
}
## OTU generate
OTU_generate<-function(city,i){
  eval(parse(text=paste0('colnames(',city,'_',i,')[1:2]<-c(\'id\',\'count\')')))
  eval(parse(text=paste0('colnames(otu_',i,')[1:3]<-c(\'id\',\'otu\',\'score\')')))
  eval(parse(text=paste0('otu_',i,'$Phyla<-NA')))
  eval(parse(text=paste0('otu_',i,'$Class<-NA')))
  eval(parse(text=paste0('otu_',i,'$Order<-NA')))
  eval(parse(text=paste0('otu_',i,'$Family<-NA')))
  eval(parse(text=paste0('otu_',i,'$Genus<-NA')))
  eval(parse(text=paste0('otu_',i,'$Species<-NA')))
  aa<-eval(parse(text=paste0('otu_',i)))
  for(j in 1:eval(parse(text=paste0('nrow(otu_',i,')')))){
    if(length(strsplit(aa[j,2],split=';')[[1]])>=2){aa$Phyla[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][2],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=3){aa$Class[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][3],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=4){aa$Order[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][4],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=5){aa$Family[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][5],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=6){aa$Genus[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2]}
    if(length(strsplit(aa[j,2],split=';')[[1]])>=7){
      aa$Species[j]<-strsplit(aa[j,2],split=';')[[1]][7]
      if(aa$Species[j]=='s__'){aa$Species[j]<-NA}
      else{aa$Species[j]<-paste0(aa$Genus[j]<-strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2],'.',strsplit(aa$Species[j],split='__')[[1]][2])}}
    if((length(strsplit(aa[j,2],split=';')[[1]])==6)|((length(strsplit(aa[j,2],split=';')[[1]])==7)&(is.na(aa$Species[j])))){
      aa$Species[j]<-paste0(strsplit(strsplit(aa[j,2],split=';')[[1]][6],split='__')[[1]][2],'.spp')
    }
    aa[which(aa$Species=='NA.spp'),'Species']<-NA
  }
  eval(parse(text=paste0('otu_',i,'<-aa')))
  eval(parse(text=paste0(city,'_',i,'<-merge(',city,'_',i,',otu_',i,',by=\'id\')')))
  aa<-get(paste0(city,'_',i))
  return(aa)
}

#############################################################################################################################
################################################ ARN: 50 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/ARN')
ARN_samples<-sprintf("%03d", 1:50)

for(i in ARN_samples){
  eval(parse(text=paste0('ARN','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in ARN_samples){
  eval(parse(text=paste0('ARN','_',j,'<-f_large_sample(','\'ARN\',\'',j,'\')')))
}

for(j in ARN_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in ARN_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in ARN_samples){
  eval(parse(text=paste0('ARN','_',i,'<-OTU_generate(','\'ARN\',\'',i,'\')')))
}

for(i in ARN_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

ARN_species<-NULL
for(i in ARN_samples){
  bb<-eval(parse(text=paste0('unique(ARN_',i,'$Species)')))
  ARN_species<-union(ARN_species,bb)
}
ARN_species<-na.omit(ARN_species)

ARN_genus<-NULL
for(i in ARN_samples){
  bb<-eval(parse(text=paste0('unique(ARN_',i,'$Genus)')))
  ARN_genus<-union(ARN_genus,bb)
}
ARN_genus<-na.omit(ARN_genus)

ARN_family<-NULL
for(i in ARN_samples){
  bb<-eval(parse(text=paste0('unique(ARN_',i,'$Family)')))
  ARN_family<-union(ARN_family,bb)
}
ARN_family<-na.omit(ARN_family)

ARN_order<-NULL
for(i in ARN_samples){
  bb<-eval(parse(text=paste0('unique(ARN_',i,'$Order)')))
  ARN_order<-union(ARN_order,bb)
}
ARN_order<-na.omit(ARN_order)

rm(i,j,bb)

#############################################################################################################################
################################################ BCN: 38 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/BCN')
BCN_samples<-sprintf("%03d", 1:38)

for(i in BCN_samples){
  eval(parse(text=paste0('BCN','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in BCN_samples){
  eval(parse(text=paste0('BCN','_',j,'<-f_large_sample(','\'BCN\',\'',j,'\')')))
}

for(j in BCN_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in BCN_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in BCN_samples){
  eval(parse(text=paste0('BCN','_',i,'<-OTU_generate(','\'BCN\',\'',i,'\')')))
}

for(i in BCN_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

BCN_species<-NULL
for(i in BCN_samples){
  bb<-eval(parse(text=paste0('unique(BCN_',i,'$Species)')))
  BCN_species<-union(BCN_species,bb)
}
BCN_species<-na.omit(BCN_species)

BCN_genus<-NULL
for(i in BCN_samples){
  bb<-eval(parse(text=paste0('unique(BCN_',i,'$Genus)')))
  BCN_genus<-union(BCN_genus,bb)
}
BCN_genus<-na.omit(BCN_genus)

BCN_family<-NULL
for(i in BCN_samples){
  bb<-eval(parse(text=paste0('unique(BCN_',i,'$Family)')))
  BCN_family<-union(BCN_family,bb)
}
BCN_family<-na.omit(BCN_family)

BCN_order<-NULL
for(i in BCN_samples){
  bb<-eval(parse(text=paste0('unique(BCN_',i,'$Order)')))
  BCN_order<-union(BCN_order,bb)
}
BCN_order<-na.omit(BCN_order)

rm(i,j,bb)

#############################################################################################################################
################################################ BER: 41 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/BER')
BER_samples<-sprintf("%03d", 1:41)

for(i in BER_samples){
  eval(parse(text=paste0('BER','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in BER_samples){
  eval(parse(text=paste0('BER','_',j,'<-f_large_sample(','\'BER\',\'',j,'\')')))
}

for(j in BER_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in BER_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in BER_samples){
  eval(parse(text=paste0('BER','_',i,'<-OTU_generate(','\'BER\',\'',i,'\')')))
}

for(i in BER_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

BER_species<-NULL
for(i in BER_samples){
  bb<-eval(parse(text=paste0('unique(BER_',i,'$Species)')))
  BER_species<-union(BER_species,bb)
}
BER_species<-na.omit(BER_species)

BER_genus<-NULL
for(i in BER_samples){
  bb<-eval(parse(text=paste0('unique(BER_',i,'$Genus)')))
  BER_genus<-union(BER_genus,bb)
}
BER_genus<-na.omit(BER_genus)

BER_family<-NULL
for(i in BER_samples){
  bb<-eval(parse(text=paste0('unique(BER_',i,'$Family)')))
  BER_family<-union(BER_family,bb)
}
BER_family<-na.omit(BER_family)

BER_order<-NULL
for(i in BER_samples){
  bb<-eval(parse(text=paste0('unique(BER_',i,'$Order)')))
  BER_order<-union(BER_order,bb)
}
BER_order<-na.omit(BER_order)

rm(i,j,bb)

#############################################################################################################################
################################################ DEN_16: 23 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/DEN_16')
DEN_16_samples<-sprintf("%03d", 1:23)

for(i in DEN_16_samples){
  eval(parse(text=paste0('DEN_16','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in DEN_16_samples){
  eval(parse(text=paste0('DEN_16','_',j,'<-f_large_sample(','\'DEN_16\',\'',j,'\')')))
}

for(j in DEN_16_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in DEN_16_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in DEN_16_samples){
  eval(parse(text=paste0('DEN_16','_',i,'<-OTU_generate(','\'DEN_16\',\'',i,'\')')))
}

for(i in DEN_16_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

DEN_16_species<-NULL
for(i in DEN_16_samples){
  bb<-eval(parse(text=paste0('unique(DEN_16_',i,'$Species)')))
  DEN_16_species<-union(DEN_16_species,bb)
}
DEN_16_species<-na.omit(DEN_16_species)

DEN_16_genus<-NULL
for(i in DEN_16_samples){
  bb<-eval(parse(text=paste0('unique(DEN_16_',i,'$Genus)')))
  DEN_16_genus<-union(DEN_16_genus,bb)
}
DEN_16_genus<-na.omit(DEN_16_genus)

DEN_16_family<-NULL
for(i in DEN_16_samples){
  bb<-eval(parse(text=paste0('unique(DEN_16_',i,'$Family)')))
  DEN_16_family<-union(DEN_16_family,bb)
}
DEN_16_family<-na.omit(DEN_16_family)

DEN_16_order<-NULL
for(i in DEN_16_samples){
  bb<-eval(parse(text=paste0('unique(DEN_16_',i,'$Order)')))
  DEN_16_order<-union(DEN_16_order,bb)
}
DEN_16_order<-na.omit(DEN_16_order)

rm(i,j,bb)

#############################################################################################################################
################################################ DEN_17: 22 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/DEN_17')
DEN_17_samples<-sprintf("%03d", 1:22)

for(i in DEN_17_samples){
  eval(parse(text=paste0('DEN_17','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in DEN_17_samples){
  eval(parse(text=paste0('DEN_17','_',j,'<-f_large_sample(','\'DEN_17\',\'',j,'\')')))
}

for(j in DEN_17_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in DEN_17_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in DEN_17_samples){
  eval(parse(text=paste0('DEN_17','_',i,'<-OTU_generate(','\'DEN_17\',\'',i,'\')')))
}

for(i in DEN_17_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

DEN_17_species<-NULL
for(i in DEN_17_samples){
  bb<-eval(parse(text=paste0('unique(DEN_17_',i,'$Species)')))
  DEN_17_species<-union(DEN_17_species,bb)
}
DEN_17_species<-na.omit(DEN_17_species)

DEN_17_genus<-NULL
for(i in DEN_17_samples){
  bb<-eval(parse(text=paste0('unique(DEN_17_',i,'$Genus)')))
  DEN_17_genus<-union(DEN_17_genus,bb)
}
DEN_17_genus<-na.omit(DEN_17_genus)

DEN_17_family<-NULL
for(i in DEN_17_samples){
  bb<-eval(parse(text=paste0('unique(DEN_17_',i,'$Family)')))
  DEN_17_family<-union(DEN_17_family,bb)
}
DEN_17_family<-na.omit(DEN_17_family)

DEN_17_order<-NULL
for(i in DEN_17_samples){
  bb<-eval(parse(text=paste0('unique(DEN_17_',i,'$Order)')))
  DEN_17_order<-union(DEN_17_order,bb)
}
DEN_17_order<-na.omit(DEN_17_order)

rm(i,j,bb)

#############################################################################################################################
################################################ DOH_16: 50 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/DOH_16')
DOH_16_samples<-sprintf("%03d", 1:48)

for(i in DOH_16_samples){
  eval(parse(text=paste0('DOH_16','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in DOH_16_samples){
  eval(parse(text=paste0('DOH_16','_',j,'<-f_large_sample(','\'DOH_16\',\'',j,'\')')))
}

for(j in DOH_16_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in DOH_16_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in DOH_16_samples){
  eval(parse(text=paste0('DOH_16','_',i,'<-OTU_generate(','\'DOH_16\',\'',i,'\')')))
}

for(i in DOH_16_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

DOH_16_species<-NULL
for(i in DOH_16_samples){
  bb<-eval(parse(text=paste0('unique(DOH_16_',i,'$Species)')))
  DOH_16_species<-union(DOH_16_species,bb)
}
DOH_16_species<-na.omit(DOH_16_species)

DOH_16_genus<-NULL
for(i in DOH_16_samples){
  bb<-eval(parse(text=paste0('unique(DOH_16_',i,'$Genus)')))
  DOH_16_genus<-union(DOH_16_genus,bb)
}
DOH_16_genus<-na.omit(DOH_16_genus)

DOH_16_family<-NULL
for(i in DOH_16_samples){
  bb<-eval(parse(text=paste0('unique(DOH_16_',i,'$Family)')))
  DOH_16_family<-union(DOH_16_family,bb)
}
DOH_16_family<-na.omit(DOH_16_family)

DOH_16_order<-NULL
for(i in DOH_16_samples){
  bb<-eval(parse(text=paste0('unique(DOH_16_',i,'$Order)')))
  DOH_16_order<-union(DOH_16_order,bb)
}
DOH_16_order<-na.omit(DOH_16_order)

rm(i,j,bb)

#############################################################################################################################
################################################ DOH_17: 15 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/DOH_17')
DOH_17_samples<-sprintf("%03d", 1:15)

for(i in DOH_17_samples){
  eval(parse(text=paste0('DOH_17','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in DOH_17_samples){
  eval(parse(text=paste0('DOH_17','_',j,'<-f_large_sample(','\'DOH_17\',\'',j,'\')')))
}

for(j in DOH_17_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in DOH_17_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in DOH_17_samples){
  eval(parse(text=paste0('DOH_17','_',i,'<-OTU_generate(','\'DOH_17\',\'',i,'\')')))
}

for(i in DOH_17_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

DOH_17_species<-NULL
for(i in DOH_17_samples){
  bb<-eval(parse(text=paste0('unique(DOH_17_',i,'$Species)')))
  DOH_17_species<-union(DOH_17_species,bb)
}
DOH_17_species<-na.omit(DOH_17_species)

DOH_17_genus<-NULL
for(i in DOH_17_samples){
  bb<-eval(parse(text=paste0('unique(DOH_17_',i,'$Genus)')))
  DOH_17_genus<-union(DOH_17_genus,bb)
}
DOH_17_genus<-na.omit(DOH_17_genus)

DOH_17_family<-NULL
for(i in DOH_17_samples){
  bb<-eval(parse(text=paste0('unique(DOH_17_',i,'$Family)')))
  DOH_17_family<-union(DOH_17_family,bb)
}
DOH_17_family<-na.omit(DOH_17_family)

DOH_17_order<-NULL
for(i in DOH_17_samples){
  bb<-eval(parse(text=paste0('unique(DOH_17_',i,'$Order)')))
  DOH_17_order<-union(DOH_17_order,bb)
}
DOH_17_order<-na.omit(DOH_17_order)

rm(i,j,bb)

#############################################################################################################################
################################################ FAI: 48 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/FAI')
FAI_samples<-sprintf("%03d", 1:48)

for(i in FAI_samples){
  eval(parse(text=paste0('FAI','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in FAI_samples){
  eval(parse(text=paste0('FAI','_',j,'<-f_large_sample(','\'FAI\',\'',j,'\')')))
}

for(j in FAI_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in FAI_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in FAI_samples){
  eval(parse(text=paste0('FAI','_',i,'<-OTU_generate(','\'FAI\',\'',i,'\')')))
}

for(i in FAI_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

FAI_species<-NULL
for(i in FAI_samples){
  bb<-eval(parse(text=paste0('unique(FAI_',i,'$Species)')))
  FAI_species<-union(FAI_species,bb)
}
FAI_species<-na.omit(FAI_species)

FAI_genus<-NULL
for(i in FAI_samples){
  bb<-eval(parse(text=paste0('unique(FAI_',i,'$Genus)')))
  FAI_genus<-union(FAI_genus,bb)
}
FAI_genus<-na.omit(FAI_genus)

FAI_family<-NULL
for(i in FAI_samples){
  bb<-eval(parse(text=paste0('unique(FAI_',i,'$Family)')))
  FAI_family<-union(FAI_family,bb)
}
FAI_family<-na.omit(FAI_family)

FAI_order<-NULL
for(i in FAI_samples){
  bb<-eval(parse(text=paste0('unique(FAI_',i,'$Order)')))
  FAI_order<-union(FAI_order,bb)
}
FAI_order<-na.omit(FAI_order)

rm(i,j,bb)

#############################################################################################################################
################################################ HKG: 49 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/HKG')
HKG_samples<-sprintf("%03d", 1:49)

for(i in HKG_samples){
  eval(parse(text=paste0('HKG','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in HKG_samples){
  eval(parse(text=paste0('HKG','_',j,'<-f_large_sample(','\'HKG\',\'',j,'\')')))
}

for(j in HKG_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in HKG_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in HKG_samples){
  eval(parse(text=paste0('HKG','_',i,'<-OTU_generate(','\'HKG\',\'',i,'\')')))
}

for(i in HKG_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

HKG_species<-NULL
for(i in HKG_samples){
  bb<-eval(parse(text=paste0('unique(HKG_',i,'$Species)')))
  HKG_species<-union(HKG_species,bb)
}
HKG_species<-na.omit(HKG_species)

HKG_genus<-NULL
for(i in HKG_samples){
  bb<-eval(parse(text=paste0('unique(HKG_',i,'$Genus)')))
  HKG_genus<-union(HKG_genus,bb)
}
HKG_genus<-na.omit(HKG_genus)

HKG_family<-NULL
for(i in HKG_samples){
  bb<-eval(parse(text=paste0('unique(HKG_',i,'$Family)')))
  HKG_family<-union(HKG_family,bb)
}
HKG_family<-na.omit(HKG_family)

HKG_order<-NULL
for(i in HKG_samples){
  bb<-eval(parse(text=paste0('unique(HKG_',i,'$Order)')))
  HKG_order<-union(HKG_order,bb)
}
HKG_order<-na.omit(HKG_order)

rm(i,j,bb)

#############################################################################################################################
################################################ ICN: 50 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/ICN')
ICN_samples<-sprintf("%03d", 1:50)

for(i in ICN_samples){
  eval(parse(text=paste0('ICN','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in ICN_samples){
  eval(parse(text=paste0('ICN','_',j,'<-f_large_sample(','\'ICN\',\'',j,'\')')))
}

for(j in ICN_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in ICN_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in ICN_samples){
  eval(parse(text=paste0('ICN','_',i,'<-OTU_generate(','\'ICN\',\'',i,'\')')))
}

for(i in ICN_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

ICN_species<-NULL
for(i in ICN_samples){
  bb<-eval(parse(text=paste0('unique(ICN_',i,'$Species)')))
  ICN_species<-union(ICN_species,bb)
}
ICN_species<-na.omit(ICN_species)

ICN_genus<-NULL
for(i in ICN_samples){
  bb<-eval(parse(text=paste0('unique(ICN_',i,'$Genus)')))
  ICN_genus<-union(ICN_genus,bb)
}
ICN_genus<-na.omit(ICN_genus)

ICN_family<-NULL
for(i in ICN_samples){
  bb<-eval(parse(text=paste0('unique(ICN_',i,'$Family)')))
  ICN_family<-union(ICN_family,bb)
}
ICN_family<-na.omit(ICN_family)

ICN_order<-NULL
for(i in ICN_samples){
  bb<-eval(parse(text=paste0('unique(ICN_',i,'$Order)')))
  ICN_order<-union(ICN_order,bb)
}
ICN_order<-na.omit(ICN_order)

rm(i,j,bb)

#############################################################################################################################
################################################ IEV: 49 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/IEV')
IEV_samples<-sprintf("%03d", 1:49)

for(i in IEV_samples){
  eval(parse(text=paste0('IEV','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in IEV_samples){
  eval(parse(text=paste0('IEV','_',j,'<-f_large_sample(','\'IEV\',\'',j,'\')')))
}

for(j in IEV_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in IEV_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in IEV_samples){
  eval(parse(text=paste0('IEV','_',i,'<-OTU_generate(','\'IEV\',\'',i,'\')')))
}

for(i in IEV_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

IEV_species<-NULL
for(i in IEV_samples){
  bb<-eval(parse(text=paste0('unique(IEV_',i,'$Species)')))
  IEV_species<-union(IEV_species,bb)
}
IEV_species<-na.omit(IEV_species)

IEV_genus<-NULL
for(i in IEV_samples){
  bb<-eval(parse(text=paste0('unique(IEV_',i,'$Genus)')))
  IEV_genus<-union(IEV_genus,bb)
}
IEV_genus<-na.omit(IEV_genus)

IEV_family<-NULL
for(i in IEV_samples){
  bb<-eval(parse(text=paste0('unique(IEV_',i,'$Family)')))
  IEV_family<-union(IEV_family,bb)
}
IEV_family<-na.omit(IEV_family)

IEV_order<-NULL
for(i in IEV_samples){
  bb<-eval(parse(text=paste0('unique(IEV_',i,'$Order)')))
  IEV_order<-union(IEV_order,bb)
}
IEV_order<-na.omit(IEV_order)

rm(i,j,bb)

#############################################################################################################################
################################################ ILR_16: 47 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/ILR_16')
ILR_16_samples<-sprintf("%03d", 1:47)

for(i in ILR_16_samples){
  eval(parse(text=paste0('ILR_16','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in ILR_16_samples){
  eval(parse(text=paste0('ILR_16','_',j,'<-f_large_sample(','\'ILR_16\',\'',j,'\')')))
}

for(j in ILR_16_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in ILR_16_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in ILR_16_samples){
  eval(parse(text=paste0('ILR_16','_',i,'<-OTU_generate(','\'ILR_16\',\'',i,'\')')))
}

for(i in ILR_16_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

ILR_16_species<-NULL
for(i in ILR_16_samples){
  bb<-eval(parse(text=paste0('unique(ILR_16_',i,'$Species)')))
  ILR_16_species<-union(ILR_16_species,bb)
}
ILR_16_species<-na.omit(ILR_16_species)

ILR_16_genus<-NULL
for(i in ILR_16_samples){
  bb<-eval(parse(text=paste0('unique(ILR_16_',i,'$Genus)')))
  ILR_16_genus<-union(ILR_16_genus,bb)
}
ILR_16_genus<-na.omit(ILR_16_genus)

ILR_16_family<-NULL
for(i in ILR_16_samples){
  bb<-eval(parse(text=paste0('unique(ILR_16_',i,'$Family)')))
  ILR_16_family<-union(ILR_16_family,bb)
}
ILR_16_family<-na.omit(ILR_16_family)

ILR_16_order<-NULL
for(i in ILR_16_samples){
  bb<-eval(parse(text=paste0('unique(ILR_16_',i,'$Order)')))
  ILR_16_order<-union(ILR_16_order,bb)
}
ILR_16_order<-na.omit(ILR_16_order)

rm(i,j,bb)

#############################################################################################################################
################################################ ILR_17: 50 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/ILR_17')
ILR_17_samples<-sprintf("%03d", 1:50)

for(i in ILR_17_samples){
  eval(parse(text=paste0('ILR_17','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in ILR_17_samples){
  eval(parse(text=paste0('ILR_17','_',j,'<-f_large_sample(','\'ILR_17\',\'',j,'\')')))
}

for(j in ILR_17_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in ILR_17_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in ILR_17_samples){
  eval(parse(text=paste0('ILR_17','_',i,'<-OTU_generate(','\'ILR_17\',\'',i,'\')')))
}

for(i in ILR_17_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

ILR_17_species<-NULL
for(i in ILR_17_samples){
  bb<-eval(parse(text=paste0('unique(ILR_17_',i,'$Species)')))
  ILR_17_species<-union(ILR_17_species,bb)
}
ILR_17_species<-na.omit(ILR_17_species)

ILR_17_genus<-NULL
for(i in ILR_17_samples){
  bb<-eval(parse(text=paste0('unique(ILR_17_',i,'$Genus)')))
  ILR_17_genus<-union(ILR_17_genus,bb)
}
ILR_17_genus<-na.omit(ILR_17_genus)

ILR_17_family<-NULL
for(i in ILR_17_samples){
  bb<-eval(parse(text=paste0('unique(ILR_17_',i,'$Family)')))
  ILR_17_family<-union(ILR_17_family,bb)
}
ILR_17_family<-na.omit(ILR_17_family)

ILR_17_order<-NULL
for(i in ILR_17_samples){
  bb<-eval(parse(text=paste0('unique(ILR_17_',i,'$Order)')))
  ILR_17_order<-union(ILR_17_order,bb)
}
ILR_17_order<-na.omit(ILR_17_order)

rm(i,j,bb)

#############################################################################################################################
################################################ KUL: 30 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/KUL')
KUL_samples<-sprintf("%03d", 1:30)

for(i in KUL_samples){
  eval(parse(text=paste0('KUL','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in KUL_samples){
  eval(parse(text=paste0('KUL','_',j,'<-f_large_sample(','\'KUL\',\'',j,'\')')))
}

for(j in KUL_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in KUL_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in KUL_samples){
  eval(parse(text=paste0('KUL','_',i,'<-OTU_generate(','\'KUL\',\'',i,'\')')))
}

for(i in KUL_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

KUL_species<-NULL
for(i in KUL_samples){
  bb<-eval(parse(text=paste0('unique(KUL_',i,'$Species)')))
  KUL_species<-union(KUL_species,bb)
}
KUL_species<-na.omit(KUL_species)

KUL_genus<-NULL
for(i in KUL_samples){
  bb<-eval(parse(text=paste0('unique(KUL_',i,'$Genus)')))
  KUL_genus<-union(KUL_genus,bb)
}
KUL_genus<-na.omit(KUL_genus)

KUL_family<-NULL
for(i in KUL_samples){
  bb<-eval(parse(text=paste0('unique(KUL_',i,'$Family)')))
  KUL_family<-union(KUL_family,bb)
}
KUL_family<-na.omit(KUL_family)

KUL_order<-NULL
for(i in KUL_samples){
  bb<-eval(parse(text=paste0('unique(KUL_',i,'$Order)')))
  KUL_order<-union(KUL_order,bb)
}
KUL_order<-na.omit(KUL_order)

rm(i,j,bb)

#############################################################################################################################
################################################ LCY: 37 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/LCY')
LCY_samples<-sprintf("%03d", 1:37)

for(i in LCY_samples){
  eval(parse(text=paste0('LCY','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in LCY_samples){
  eval(parse(text=paste0('LCY','_',j,'<-f_large_sample(','\'LCY\',\'',j,'\')')))
}

for(j in LCY_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in LCY_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in LCY_samples){
  eval(parse(text=paste0('LCY','_',i,'<-OTU_generate(','\'LCY\',\'',i,'\')')))
}

for(i in LCY_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

LCY_species<-NULL
for(i in LCY_samples){
  bb<-eval(parse(text=paste0('unique(LCY_',i,'$Species)')))
  LCY_species<-union(LCY_species,bb)
}
LCY_species<-na.omit(LCY_species)

LCY_genus<-NULL
for(i in LCY_samples){
  bb<-eval(parse(text=paste0('unique(LCY_',i,'$Genus)')))
  LCY_genus<-union(LCY_genus,bb)
}
LCY_genus<-na.omit(LCY_genus)

LCY_family<-NULL
for(i in LCY_samples){
  bb<-eval(parse(text=paste0('unique(LCY_',i,'$Family)')))
  LCY_family<-union(LCY_family,bb)
}
LCY_family<-na.omit(LCY_family)

LCY_order<-NULL
for(i in LCY_samples){
  bb<-eval(parse(text=paste0('unique(LCY_',i,'$Order)')))
  LCY_order<-union(LCY_order,bb)
}
LCY_order<-na.omit(LCY_order)

rm(i,j,bb)

#############################################################################################################################
################################################ LIS: 19 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/LIS')
LIS_samples<-sprintf("%03d", 1:19)

for(i in LIS_samples){
  eval(parse(text=paste0('LIS','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in LIS_samples){
  eval(parse(text=paste0('LIS','_',j,'<-f_large_sample(','\'LIS\',\'',j,'\')')))
}

for(j in LIS_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in LIS_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in LIS_samples){
  eval(parse(text=paste0('LIS','_',i,'<-OTU_generate(','\'LIS\',\'',i,'\')')))
}

for(i in LIS_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

LIS_species<-NULL
for(i in LIS_samples){
  bb<-eval(parse(text=paste0('unique(LIS_',i,'$Species)')))
  LIS_species<-union(LIS_species,bb)
}
LIS_species<-na.omit(LIS_species)

LIS_genus<-NULL
for(i in LIS_samples){
  bb<-eval(parse(text=paste0('unique(LIS_',i,'$Genus)')))
  LIS_genus<-union(LIS_genus,bb)
}
LIS_genus<-na.omit(LIS_genus)

LIS_family<-NULL
for(i in LIS_samples){
  bb<-eval(parse(text=paste0('unique(LIS_',i,'$Family)')))
  LIS_family<-union(LIS_family,bb)
}
LIS_family<-na.omit(LIS_family)

LIS_order<-NULL
for(i in LIS_samples){
  bb<-eval(parse(text=paste0('unique(LIS_',i,'$Order)')))
  LIS_order<-union(LIS_order,bb)
}
LIS_order<-na.omit(LIS_order)

rm(i,j,bb)

#############################################################################################################################
################################################ NYC_16: 49 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/NYC_16')
NYC_16_samples<-sprintf("%03d", 1:49)

for(i in NYC_16_samples){
  eval(parse(text=paste0('NYC_16','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in NYC_16_samples){
  eval(parse(text=paste0('NYC_16','_',j,'<-f_large_sample(','\'NYC_16\',\'',j,'\')')))
}

for(j in NYC_16_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in NYC_16_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in NYC_16_samples){
  eval(parse(text=paste0('NYC_16','_',i,'<-OTU_generate(','\'NYC_16\',\'',i,'\')')))
}

for(i in NYC_16_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

NYC_16_species<-NULL
for(i in NYC_16_samples){
  bb<-eval(parse(text=paste0('unique(NYC_16_',i,'$Species)')))
  NYC_16_species<-union(NYC_16_species,bb)
}
NYC_16_species<-na.omit(NYC_16_species)

NYC_16_genus<-NULL
for(i in NYC_16_samples){
  bb<-eval(parse(text=paste0('unique(NYC_16_',i,'$Genus)')))
  NYC_16_genus<-union(NYC_16_genus,bb)
}
NYC_16_genus<-na.omit(NYC_16_genus)

NYC_16_family<-NULL
for(i in NYC_16_samples){
  bb<-eval(parse(text=paste0('unique(NYC_16_',i,'$Family)')))
  NYC_16_family<-union(NYC_16_family,bb)
}
NYC_16_family<-na.omit(NYC_16_family)

NYC_16_order<-NULL
for(i in NYC_16_samples){
  bb<-eval(parse(text=paste0('unique(NYC_16_',i,'$Order)')))
  NYC_16_order<-union(NYC_16_order,bb)
}
NYC_16_order<-na.omit(NYC_16_order)

rm(i,j,bb)

#############################################################################################################################
################################################ NYC_17: 50 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/NYC_17')
NYC_17_samples<-sprintf("%03d", 1:50)

for(i in NYC_17_samples){
  eval(parse(text=paste0('NYC_17','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in NYC_17_samples){
  eval(parse(text=paste0('NYC_17','_',j,'<-f_large_sample(','\'NYC_17\',\'',j,'\')')))
}

for(j in NYC_17_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in NYC_17_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in NYC_17_samples){
  eval(parse(text=paste0('NYC_17','_',i,'<-OTU_generate(','\'NYC_17\',\'',i,'\')')))
}

for(i in NYC_17_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

NYC_17_species<-NULL
for(i in NYC_17_samples){
  bb<-eval(parse(text=paste0('unique(NYC_17_',i,'$Species)')))
  NYC_17_species<-union(NYC_17_species,bb)
}
NYC_17_species<-na.omit(NYC_17_species)

NYC_17_genus<-NULL
for(i in NYC_17_samples){
  bb<-eval(parse(text=paste0('unique(NYC_17_',i,'$Genus)')))
  NYC_17_genus<-union(NYC_17_genus,bb)
}
NYC_17_genus<-na.omit(NYC_17_genus)

NYC_17_family<-NULL
for(i in NYC_17_samples){
  bb<-eval(parse(text=paste0('unique(NYC_17_',i,'$Family)')))
  NYC_17_family<-union(NYC_17_family,bb)
}
NYC_17_family<-na.omit(NYC_17_family)

NYC_17_order<-NULL
for(i in NYC_17_samples){
  bb<-eval(parse(text=paste0('unique(NYC_17_',i,'$Order)')))
  NYC_17_order<-union(NYC_17_order,bb)
}
NYC_17_order<-na.omit(NYC_17_order)

rm(i,j,bb)

#############################################################################################################################
################################################ OFF: 26 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/OFF')
OFF_samples<-sprintf("%03d", 1:26)

for(i in OFF_samples){
  eval(parse(text=paste0('OFF','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in OFF_samples){
  eval(parse(text=paste0('OFF','_',j,'<-f_large_sample(','\'OFF\',\'',j,'\')')))
}

for(j in OFF_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in OFF_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in OFF_samples){
  eval(parse(text=paste0('OFF','_',i,'<-OTU_generate(','\'OFF\',\'',i,'\')')))
}

for(i in OFF_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

OFF_species<-NULL
for(i in OFF_samples){
  bb<-eval(parse(text=paste0('unique(OFF_',i,'$Species)')))
  OFF_species<-union(OFF_species,bb)
}
OFF_species<-na.omit(OFF_species)

OFF_genus<-NULL
for(i in OFF_samples){
  bb<-eval(parse(text=paste0('unique(OFF_',i,'$Genus)')))
  OFF_genus<-union(OFF_genus,bb)
}
OFF_genus<-na.omit(OFF_genus)

OFF_family<-NULL
for(i in OFF_samples){
  bb<-eval(parse(text=paste0('unique(OFF_',i,'$Family)')))
  OFF_family<-union(OFF_family,bb)
}
OFF_family<-na.omit(OFF_family)

OFF_order<-NULL
for(i in OFF_samples){
  bb<-eval(parse(text=paste0('unique(OFF_',i,'$Order)')))
  OFF_order<-union(OFF_order,bb)
}
OFF_order<-na.omit(OFF_order)

rm(i,j,bb)

#############################################################################################################################
################################################ SAO: 29 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/SAO')
SAO_samples<-sprintf("%03d", 1:29)

for(i in SAO_samples){
  eval(parse(text=paste0('SAO','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in SAO_samples){
  eval(parse(text=paste0('SAO','_',j,'<-f_large_sample(','\'SAO\',\'',j,'\')')))
}

for(j in SAO_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in SAO_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in SAO_samples){
  eval(parse(text=paste0('SAO','_',i,'<-OTU_generate(','\'SAO\',\'',i,'\')')))
}

for(i in SAO_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

SAO_species<-NULL
for(i in SAO_samples){
  bb<-eval(parse(text=paste0('unique(SAO_',i,'$Species)')))
  SAO_species<-union(SAO_species,bb)
}
SAO_species<-na.omit(SAO_species)

SAO_genus<-NULL
for(i in SAO_samples){
  bb<-eval(parse(text=paste0('unique(SAO_',i,'$Genus)')))
  SAO_genus<-union(SAO_genus,bb)
}
SAO_genus<-na.omit(SAO_genus)

SAO_family<-NULL
for(i in SAO_samples){
  bb<-eval(parse(text=paste0('unique(SAO_',i,'$Family)')))
  SAO_family<-union(SAO_family,bb)
}
SAO_family<-na.omit(SAO_family)

SAO_order<-NULL
for(i in SAO_samples){
  bb<-eval(parse(text=paste0('unique(SAO_',i,'$Order)')))
  SAO_order<-union(SAO_order,bb)
}
SAO_order<-na.omit(SAO_order)

rm(i,j,bb)

#############################################################################################################################
################################################ SCL: 26 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/SCL')
SCL_samples<-sprintf("%03d", 1:26)

for(i in SCL_samples){
  eval(parse(text=paste0('SCL','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in SCL_samples){
  eval(parse(text=paste0('SCL','_',j,'<-f_large_sample(','\'SCL\',\'',j,'\')')))
}

for(j in SCL_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in SCL_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in SCL_samples){
  eval(parse(text=paste0('SCL','_',i,'<-OTU_generate(','\'SCL\',\'',i,'\')')))
}

for(i in SCL_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

SCL_species<-NULL
for(i in SCL_samples){
  bb<-eval(parse(text=paste0('unique(SCL_',i,'$Species)')))
  SCL_species<-union(SCL_species,bb)
}
SCL_species<-na.omit(SCL_species)

SCL_genus<-NULL
for(i in SCL_samples){
  bb<-eval(parse(text=paste0('unique(SCL_',i,'$Genus)')))
  SCL_genus<-union(SCL_genus,bb)
}
SCL_genus<-na.omit(SCL_genus)

SCL_family<-NULL
for(i in SCL_samples){
  bb<-eval(parse(text=paste0('unique(SCL_',i,'$Family)')))
  SCL_family<-union(SCL_family,bb)
}
SCL_family<-na.omit(SCL_family)

SCL_order<-NULL
for(i in SCL_samples){
  bb<-eval(parse(text=paste0('unique(SCL_',i,'$Order)')))
  SCL_order<-union(SCL_order,bb)
}
SCL_order<-na.omit(SCL_order)

rm(i,j,bb)

#############################################################################################################################
################################################ SDJ: 32 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/SDJ')
SDJ_samples<-sprintf("%03d", 1:32)

for(i in SDJ_samples){
  eval(parse(text=paste0('SDJ','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in SDJ_samples){
  eval(parse(text=paste0('SDJ','_',j,'<-f_large_sample(','\'SDJ\',\'',j,'\')')))
}

for(j in SDJ_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in SDJ_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in SDJ_samples){
  eval(parse(text=paste0('SDJ','_',i,'<-OTU_generate(','\'SDJ\',\'',i,'\')')))
}

for(i in SDJ_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

SDJ_species<-NULL
for(i in SDJ_samples){
  bb<-eval(parse(text=paste0('unique(SDJ_',i,'$Species)')))
  SDJ_species<-union(SDJ_species,bb)
}
SDJ_species<-na.omit(SDJ_species)

SDJ_genus<-NULL
for(i in SDJ_samples){
  bb<-eval(parse(text=paste0('unique(SDJ_',i,'$Genus)')))
  SDJ_genus<-union(SDJ_genus,bb)
}
SDJ_genus<-na.omit(SDJ_genus)

SDJ_family<-NULL
for(i in SDJ_samples){
  bb<-eval(parse(text=paste0('unique(SDJ_',i,'$Family)')))
  SDJ_family<-union(SDJ_family,bb)
}
SDJ_family<-na.omit(SDJ_family)

SDJ_order<-NULL
for(i in SDJ_samples){
  bb<-eval(parse(text=paste0('unique(SDJ_',i,'$Order)')))
  SDJ_order<-union(SDJ_order,bb)
}
SDJ_order<-na.omit(SDJ_order)

rm(i,j,bb)

#############################################################################################################################
################################################ SFO: 27 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/SFO')
SFO_samples<-sprintf("%03d", 1:29)
SFO_samples<-setdiff(SFO_samples,c('010','025'))

for(i in SFO_samples){
  eval(parse(text=paste0('SFO','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in SFO_samples){
  eval(parse(text=paste0('SFO','_',j,'<-f_large_sample(','\'SFO\',\'',j,'\')')))
}

for(j in SFO_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in SFO_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in SFO_samples){
  eval(parse(text=paste0('SFO','_',i,'<-OTU_generate(','\'SFO\',\'',i,'\')')))
}

for(i in SFO_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

SFO_species<-NULL
for(i in SFO_samples){
  bb<-eval(parse(text=paste0('unique(SFO_',i,'$Species)')))
  SFO_species<-union(SFO_species,bb)
}
SFO_species<-na.omit(SFO_species)

SFO_genus<-NULL
for(i in SFO_samples){
  bb<-eval(parse(text=paste0('unique(SFO_',i,'$Genus)')))
  SFO_genus<-union(SFO_genus,bb)
}
SFO_genus<-na.omit(SFO_genus)

SFO_family<-NULL
for(i in SFO_samples){
  bb<-eval(parse(text=paste0('unique(SFO_',i,'$Family)')))
  SFO_family<-union(SFO_family,bb)
}
SFO_family<-na.omit(SFO_family)

SFO_order<-NULL
for(i in SFO_samples){
  bb<-eval(parse(text=paste0('unique(SFO_',i,'$Order)')))
  SFO_order<-union(SFO_order,bb)
}
SFO_order<-na.omit(SFO_order)

rm(i,j,bb)

#############################################################################################################################
################################################ SGP: 48 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/SGP')
SGP_samples<-sprintf("%03d", 1:48)

for(i in SGP_samples){
  eval(parse(text=paste0('SGP','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in SGP_samples){
  eval(parse(text=paste0('SGP','_',j,'<-f_large_sample(','\'SGP\',\'',j,'\')')))
}

for(j in SGP_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in SGP_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in SGP_samples){
  eval(parse(text=paste0('SGP','_',i,'<-OTU_generate(','\'SGP\',\'',i,'\')')))
}

for(i in SGP_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

SGP_species<-NULL
for(i in SGP_samples){
  bb<-eval(parse(text=paste0('unique(SGP_',i,'$Species)')))
  SGP_species<-union(SGP_species,bb)
}
SGP_species<-na.omit(SGP_species)

SGP_genus<-NULL
for(i in SGP_samples){
  bb<-eval(parse(text=paste0('unique(SGP_',i,'$Genus)')))
  SGP_genus<-union(SGP_genus,bb)
}
SGP_genus<-na.omit(SGP_genus)

SGP_family<-NULL
for(i in SGP_samples){
  bb<-eval(parse(text=paste0('unique(SGP_',i,'$Family)')))
  SGP_family<-union(SGP_family,bb)
}
SGP_family<-na.omit(SGP_family)

SGP_order<-NULL
for(i in SGP_samples){
  bb<-eval(parse(text=paste0('unique(SGP_',i,'$Order)')))
  SGP_order<-union(SGP_order,bb)
}
SGP_order<-na.omit(SGP_order)

rm(i,j,bb)

#############################################################################################################################
################################################ TPE: 50 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/TPE')
TPE_samples<-sprintf("%03d", 1:50)

for(i in TPE_samples){
  eval(parse(text=paste0('TPE','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in TPE_samples){
  eval(parse(text=paste0('TPE','_',j,'<-f_large_sample(','\'TPE\',\'',j,'\')')))
}

for(j in TPE_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in TPE_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in TPE_samples){
  eval(parse(text=paste0('TPE','_',i,'<-OTU_generate(','\'TPE\',\'',i,'\')')))
}

for(i in TPE_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

TPE_species<-NULL
for(i in TPE_samples){
  bb<-eval(parse(text=paste0('unique(TPE_',i,'$Species)')))
  TPE_species<-union(TPE_species,bb)
}
TPE_species<-na.omit(TPE_species)

TPE_genus<-NULL
for(i in TPE_samples){
  bb<-eval(parse(text=paste0('unique(TPE_',i,'$Genus)')))
  TPE_genus<-union(TPE_genus,bb)
}
TPE_genus<-na.omit(TPE_genus)

TPE_family<-NULL
for(i in TPE_samples){
  bb<-eval(parse(text=paste0('unique(TPE_',i,'$Family)')))
  TPE_family<-union(TPE_family,bb)
}
TPE_family<-na.omit(TPE_family)

TPE_order<-NULL
for(i in TPE_samples){
  bb<-eval(parse(text=paste0('unique(TPE_',i,'$Order)')))
  TPE_order<-union(TPE_order,bb)
}
TPE_order<-na.omit(TPE_order)

rm(i,j,bb)

#############################################################################################################################
################################################ TYO_16: 25 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/TYO_16')
TYO_16_samples<-sprintf("%03d", 1:25)

for(i in TYO_16_samples){
  eval(parse(text=paste0('TYO_16','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in TYO_16_samples){
  eval(parse(text=paste0('TYO_16','_',j,'<-f_large_sample(','\'TYO_16\',\'',j,'\')')))
}

for(j in TYO_16_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in TYO_16_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in TYO_16_samples){
  eval(parse(text=paste0('TYO_16','_',i,'<-OTU_generate(','\'TYO_16\',\'',i,'\')')))
}

for(i in TYO_16_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

TYO_16_species<-NULL
for(i in TYO_16_samples){
  bb<-eval(parse(text=paste0('unique(TYO_16_',i,'$Species)')))
  TYO_16_species<-union(TYO_16_species,bb)
}
TYO_16_species<-na.omit(TYO_16_species)

TYO_16_genus<-NULL
for(i in TYO_16_samples){
  bb<-eval(parse(text=paste0('unique(TYO_16_',i,'$Genus)')))
  TYO_16_genus<-union(TYO_16_genus,bb)
}
TYO_16_genus<-na.omit(TYO_16_genus)

TYO_16_family<-NULL
for(i in TYO_16_samples){
  bb<-eval(parse(text=paste0('unique(TYO_16_',i,'$Family)')))
  TYO_16_family<-union(TYO_16_family,bb)
}
TYO_16_family<-na.omit(TYO_16_family)

TYO_16_order<-NULL
for(i in TYO_16_samples){
  bb<-eval(parse(text=paste0('unique(TYO_16_',i,'$Order)')))
  TYO_16_order<-union(TYO_16_order,bb)
}
TYO_16_order<-na.omit(TYO_16_order)

rm(i,j,bb)

#############################################################################################################################
################################################ TYO_17: 49 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/TYO_17')
TYO_17_samples<-sprintf("%03d", 1:50)
TYO_17_samples<-setdiff(TYO_17_samples,'030')

for(i in TYO_17_samples){
  eval(parse(text=paste0('TYO_17','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in TYO_17_samples){
  eval(parse(text=paste0('TYO_17','_',j,'<-f_large_sample(','\'TYO_17\',\'',j,'\')')))
}

for(j in TYO_17_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in TYO_17_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in TYO_17_samples){
  eval(parse(text=paste0('TYO_17','_',i,'<-OTU_generate(','\'TYO_17\',\'',i,'\')')))
}

for(i in TYO_17_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

TYO_17_species<-NULL
for(i in TYO_17_samples){
  bb<-eval(parse(text=paste0('unique(TYO_17_',i,'$Species)')))
  TYO_17_species<-union(TYO_17_species,bb)
}
TYO_17_species<-na.omit(TYO_17_species)

TYO_17_genus<-NULL
for(i in TYO_17_samples){
  bb<-eval(parse(text=paste0('unique(TYO_17_',i,'$Genus)')))
  TYO_17_genus<-union(TYO_17_genus,bb)
}
TYO_17_genus<-na.omit(TYO_17_genus)

TYO_17_family<-NULL
for(i in TYO_17_samples){
  bb<-eval(parse(text=paste0('unique(TYO_17_',i,'$Family)')))
  TYO_17_family<-union(TYO_17_family,bb)
}
TYO_17_family<-na.omit(TYO_17_family)

TYO_17_order<-NULL
for(i in TYO_17_samples){
  bb<-eval(parse(text=paste0('unique(TYO_17_',i,'$Order)')))
  TYO_17_order<-union(TYO_17_order,bb)
}
TYO_17_order<-na.omit(TYO_17_order)

rm(i,j,bb)

#############################################################################################################################
################################################ ZRH: 33 samples#############################################################
setwd('E:/University of Florida/My research/CAMDA 2020/data/ZRH')
ZRH_samples<-sprintf("%03d", 1:33)

for(i in ZRH_samples){
  eval(parse(text=paste0('ZRH','_',i,'<-read.delim(\'',i,'/wgs_from_biom_final.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
  eval(parse(text=paste0('otu_',i,'<-read.delim(\'',i,'/rep_set_aligned_tax_assignments.txt\', header=FALSE, comment.char=\'#\',stringsAsFactors=FALSE)')))
}

for(j in ZRH_samples){
  eval(parse(text=paste0('ZRH','_',j,'<-f_large_sample(','\'ZRH\',\'',j,'\')')))
}

for(j in ZRH_samples){
  otu_info<-eval(parse(text=paste0('otu_',j)))
  if(ncol(otu_info)==4){
    bb<-otu_info[!duplicated(otu_info$V2),-1]
    eval(parse(text=paste0('otu_',j,'<-bb')))
    rm(bb)
  }
  rm(otu_info)
}

for(i in ZRH_samples){
  eval(parse(text=paste0('otu_',i,'[,2]<-gsub(\'[][]\',\'\',otu_',i,'[,2])')))
}

for(i in ZRH_samples){
  eval(parse(text=paste0('ZRH','_',i,'<-OTU_generate(','\'ZRH\',\'',i,'\')')))
}

for(i in ZRH_samples){
  eval(parse(text=paste0('rm(otu_',i,')')))
}

ZRH_species<-NULL
for(i in ZRH_samples){
  bb<-eval(parse(text=paste0('unique(ZRH_',i,'$Species)')))
  ZRH_species<-union(ZRH_species,bb)
}
ZRH_species<-na.omit(ZRH_species)

ZRH_genus<-NULL
for(i in ZRH_samples){
  bb<-eval(parse(text=paste0('unique(ZRH_',i,'$Genus)')))
  ZRH_genus<-union(ZRH_genus,bb)
}
ZRH_genus<-na.omit(ZRH_genus)

ZRH_family<-NULL
for(i in ZRH_samples){
  bb<-eval(parse(text=paste0('unique(ZRH_',i,'$Family)')))
  ZRH_family<-union(ZRH_family,bb)
}
ZRH_family<-na.omit(ZRH_family)

ZRH_order<-NULL
for(i in ZRH_samples){
  bb<-eval(parse(text=paste0('unique(ZRH_',i,'$Order)')))
  ZRH_order<-union(ZRH_order,bb)
}
ZRH_order<-na.omit(ZRH_order)

rm(i,j,bb)

###### OUTPUT ########
###### OUTPUT ########
###### OUTPUT ########
###### OUTPUT ########
setwd('E:/University of Florida/My research/CAMDA 2020/data')
save.image('Main dataset 1 Generate OTU tables for all samples.RData')
all_variables<-ls()
write.table(all_variables,'all_variables.txt')


