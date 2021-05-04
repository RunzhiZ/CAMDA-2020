setwd('E:/University of Florida/My research/CAMDA 2020/data')
library(limma)
library(edgeR)
library(vegan)
library(ggplot2)
library(ape)


data_for_PCoA<-data_used_normalized_12

#### Vector of cities ####
#### Vector of cities ####
#### Vector of cities ####
PCoA_city <- data_for_PCoA[,'city']
PCoA_otu<- data.frame((data_for_PCoA[,-ncol(data_for_PCoA)]))

#### PCoA ####
#### PCoA ####
#### PCoA ####
dist <- vegdist(PCoA_otu,  method = "bray",correction = "cailliez")
PCOA <- pcoa(dist)
PCOA$values$Relative_eig

pcoa_plot<-as.data.frame(matrix(0,nrow=nrow(PCoA_otu),ncol=4))
pcoa_plot[,1]<-row.names(PCoA_otu)
pcoa_plot[,2]<-PCoA_city
pcoa_plot[,3]<-PCOA$vectors[,1]
pcoa_plot[,4]<-PCOA$vectors[,2]
colnames(pcoa_plot)<-c('sample','group','PCoA1','PCoA2')

p_plot<-ggplot(pcoa_plot,aes(PCoA1,PCoA2,group=group,color=group))+
  geom_point()+
  stat_ellipse(level = 0.95, show.legend = F)+
  labs(x = paste('PCoA1: ',round(PCOA$values$Relative_eig[1]*100,1),'%'), y = paste('PCoA2: ',round(PCOA$values$Relative_eig[2]*100,1),'%'))+
  guides(colour = guide_legend(nrow = 2, override.aes = list(size=1) )) + labs(title = paste0("PCoA - Normalized data\n","Main dataset","\nCity clusters") ) + 
  theme(panel.background=element_rect(fill=NA),
        panel.border = element_rect(fill = NA, colour="black", size = 0.7),
        panel.grid.major = element_line(size=0.35, colour = "black", linetype=3),
        panel.grid.minor = element_line(size=0.15, colour = "gray50", linetype=3),
        legend.title=element_blank(),
        legend.text.align = 0, legend.text=element_text(size=5),
        legend.direction = "horizontal", legend.position = "bottom",
        legend.key = element_rect(fill = "white", colour = "white"), legend.key.size = unit(1,"line"),
        plot.title = element_text(lineheight=1, size = 7 ,face="bold"),
        axis.text.x= element_text(colour = "black", size = 8, angle = 0, vjust = 0.5),
        axis.text.y = element_text(colour = "black", size = 8),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10, angle = 90, vjust = 0.3)) #+ xlim(-5,5) + ylim(-5,5)
print(p_plot)
tiff(paste0( 'PCoA main dataset',".tif"), height=3000, width=5000, res=600, units="px"); print(p_plot); dev.off()
