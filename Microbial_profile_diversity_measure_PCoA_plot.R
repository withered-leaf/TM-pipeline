#Construction of microbiome profile 
#Written by Wei Guo
#Data:Nov 20,2020

#load library
library(matrixStats)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(plyr)
library(ggplot2)
library(ggpubr)
library(car)
library(vegan)
library(egg)
#load microbial abundance matrix from kraken2+Bracken output.
setwd("/data1/PDAC_files/code_ocean/")
#the_taxonomic_level="genus"
myfile="PDAC_samples_count_genus"
bac.exp<-read.delim(file = myfile,row.names = 1,stringsAsFactors = F)
#relative abundance matrix
bac.exp<-as.matrix(bac.exp)
sum_matrix<-matrix(rep(colSums2(bac.exp),nrow(bac.exp)),nrow=nrow(bac.exp),byrow=T)
rela.bac.exp<-bac.exp/sum_matrix
#### load group information
Groupinfo<-read.delim("Sample_group.txt",stringsAsFactors = F)
Groupinfo$Subtype<-factor(Groupinfo$Subtype,levels = c("Classical","Hybrid","Basal"))
bac.exp<-rela.bac.exp[,Groupinfo$Samples]
#plot predominant microorganisms by group
groups_mean.exp<-apply(bac.exp,1,function(x){
  groups1.mean=mean(x[which(Groupinfo$Subtype=="Classical")])
  groups2.mean=mean(x[which(Groupinfo$Subtype=="Hybrid")])
  groups3.mean=mean(x[which(Groupinfo$Subtype=="Basal")])
  merge.mean=cbind(groups1.mean,groups2.mean,groups3.mean)
  return(merge.mean)
})
groups_mean.exp<-t(groups_mean.exp)
colnames(groups_mean.exp)=c('Classical','Hybrid','Basal')
micro.sum=apply(groups_mean.exp,1,sum)
#sort by the abundance
res1=cbind(groups_mean.exp,sum=micro.sum)
res2=res1[order(-res1[,4]),]
#pick the top candidates for examples: top 15
candidates.list<-row.names(res2[1:15,])
candidates.bac.rela.exp<-groups_mean.exp[candidates.list,]
supplyothers=data.frame(Classical=1-colSums2(candidates.bac.rela.exp)[1],Hybrid=1-colSums2(candidates.bac.rela.exp)[2],Basal=1-colSums2(candidates.bac.rela.exp)[3])
res3=rbind(candidates.bac.rela.exp,others=supplyothers)
res4<-melt(as.matrix(res3))
colnames(res4)=c("genus","group","counts")
mycol<-colorRampPalette(brewer.pal(12, "Set3"))(16)
res=ggplot(data=res4,aes(x=group,y=counts,fill=genus))+geom_bar(stat="identity",color='black')+
  scale_fill_manual(values = mycol)+ylab("relative abundance")+ggtitle("Major microorganisms profile")+
  xlab("")+theme(axis.text.x = element_text(size = 8))
ggsave(filename = "predominant_genus_barplot_eachgroup.pdf",res,width = 4.5,height = 5.5)

#plot predominant microorganisms by sample
candidates.bac.rela.data<-bac.exp[candidates.list,]
supplepthers<- 1- colSums2(candidates.bac.rela.data)
res3<-rbind(candidates.bac.rela.data,Others=supplepthers)
res4<-melt(as.matrix(res3))
colnames(res4)=c("genus","sample","counts")
res4$genus<-factor(res4$genus,levels = row.names(res3))
res4$group<-factor(Groupinfo$Subtype[match(res4$sample,Groupinfo$Samples)])
res=ggplot(res4 %>% arrange(-counts) %>% 
             mutate(sample=factor(sample,levels=sample[genus=="Pseudomonas"]))
           ,aes(x=sample,y=counts,fill=genus),split="group")+geom_bar(stat="identity",color='black')+
  scale_fill_manual(values = mycol)+ylab("relative abundance")+ggtitle("")+
  xlab("")+facet_grid(. ~ group,scales="free", space="free_x")+
  theme(axis.text.x = element_blank(),axis.title.x=element_blank(),axis.ticks.x=element_blank())
ggsave(filename = "predominant_genus_barplot_eachsample_facet.pdf",res,width = 8,height = 5.5)

##microbial diversity measurement
#calculate shannon-index for the relative abundance file. default: row=microbiobe  col=samples.
mymatrix<-t(rela.bac.exp)
shannon<-function(x){sum(x[which(x>0)]*log(x[which(x>0)])*-1)}
index<-apply(mymatrix,1,shannon)
index<-data.frame(shannon=index,group=Groupinfo$Subtype,type="Shannon")
qqPlot(lm(shannon ~ group, data = index), simulate = TRUE, main = "QQ Plot", labels = FALSE)
fit_shannon<-aov(shannon~group,data=index)
TukeyHSD(fit_shannon)
my_comparison<-list(c("Classical","Hybrid"),c("Hybrid","Basal"),c("Classical","Basal"))
p1<-ggboxplot(index,x="group",y="shannon",color="group",outlier.shape = NA)+stat_compare_means(comparisons = my_comparison,method = "t.test")+
  stat_compare_means(method = "anova")+
  scale_color_manual(values=c("deepskyblue","gray20","brown1"))+
  geom_jitter(aes(color=group),size=0.8)+
  ggtitle("Shannon")+ylab("")+xlab("")+guides(color=F)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x =element_blank(),legend.title=element_blank(),
        panel.background = element_rect(fill = "grey90",colour = NA), panel.border = element_blank(),
        panel.grid.major = element_line(colour = "white"),axis.line = element_blank())
p1
#calculate richness for the relative abundance file.
richness<-function(x){length(x[which(x>0)])}
r_index<-apply(mymatrix,1,richness)
r_index<-data.frame(richness=r_index,group=Groupinfo$Subtype,type="Richness")
fit_richness<-aov(richness~group,data=r_index)
TukeyHSD(fit_richness)
p2<-ggboxplot(r_index,x="group",y="richness",color="group",outlier.shape = NA)+stat_compare_means(comparisons = my_comparison,method = "t.test")+
  stat_compare_means(method="anova")+
  scale_color_manual(values=c("deepskyblue","gray20","brown1"))+
  geom_jitter(aes(color=group),size=0.8)+
  ggtitle("Richness")+ylab("Diversity Measure")+xlab("")+guides(color=F)+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x =element_blank(),legend.title=element_blank(),
        panel.background = element_rect(fill = "grey90",colour = NA), panel.border = element_blank(),
        panel.grid.major = element_line(colour = "white"),axis.line = element_blank())
p2

#calculate bray-curtis distance 
distance <- vegdist(mymatrix, method = 'bray')
distance<-as.matrix(distance)
calculate.b.div<-function(x){
  nums<-dim(x)[1]
  out<-c()
  for(i in 1:(nums-1)){
    tmp_out<-x[,i][(i+1):nums]
    out<-c(out,tmp_out)
  }
  return(out)
}
calculate.b.div.within.groups<-function(x,Groupinfo,group){
  index<-grep(group,Groupinfo$Subtype)
  my_mtx<-x[index,index]
  out<-calculate.b.div(my_mtx)
  return(out)
}
b.div.group1<-calculate.b.div.within.groups(distance,Groupinfo,"Classical")
b.div.group2<-calculate.b.div.within.groups(distance,Groupinfo,"Hybrid")
b.div.group3<-calculate.b.div.within.groups(distance,Groupinfo,"Basal")
b_index<-data.frame(bray=c(b.div.group1,b.div.group2,b.div.group3),group=c(rep("Classical",time=length(b.div.group1)),rep("Hybrid",time=length(b.div.group2)),rep("Basal",time=length(b.div.group3))),type="Bray-Crutis Distance")
b_index$group<-factor(b_index$group,levels = c("Classical","Hybrid","Basal"))
fit_bray<-aov(bray~group,data=b_index)
TukeyHSD(fit_bray)
p3<-ggboxplot(b_index,x="group",y="bray",color="group",outlier.shape = NA)+stat_compare_means(comparisons = my_comparison,method = "t.test")+
  stat_compare_means(method="anova")+
  scale_color_manual(values=c("deepskyblue","gray20","brown1"))+
  #geom_jitter(aes(color=group),size=0.8)+
  ggtitle("Bray-Crutis Distance")+ylab("")+xlab("")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.x =element_blank(),legend.title=element_blank(),
        panel.background = element_rect(fill = "grey90",colour = NA), panel.border = element_blank(),
        panel.grid.major = element_line(colour = "white"),axis.line = element_blank(),
        legend.position = "right")
p3

egg::ggarrange(p2,p1,p3,ncol=3)

###plot PCoA between three groups
pcoa <- cmdscale(distance, k = (nrow(mymatrix) - 1), eig = TRUE)
#Axis interpretation (first two axes)
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
sample_site <- data.frame({pcoa$point})[1:2]
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')
#Add grouping information for sample point coordinates
sample_site <- merge(sample_site, Groupinfo, by.x = 'names',by.y = "Samples", all.x = TRUE)
colnames(sample_site)[4]<-"group"
sample_site$group <- factor(sample_site$group, levels = c("Classical","Hybrid","Basal"))
#
#ellipse  
pcoa_p1<-ggplot(sample_site,aes(PCoA1,PCoA2, color = group))+
  geom_point(aes(color = group), size = 1.5, alpha = 0.7)+
  stat_ellipse(aes( PCoA1,PCoA2,fill=group),geom="polygon",level=0.9,alpha=0.1)+
  scale_fill_manual(values=c("deepskyblue","gray20","brown1"))+
  scale_color_manual(values = c("deepskyblue","gray20","brown1"))+
  labs(x = paste('PCoA axis1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA axis2: ', round(100 * pcoa_eig[2], 2), '%'))+
  theme(legend.title = element_blank())
pcoa_p1 

# PERMANOVA analyse
# 
adonis_result_dis = adonis(distance~group, sample_site, permutations = 999)
adonis_result_dis
# paired comparison
group_name = unique(sample_site$group)
result = data.frame()
for (i in 1:(length(group_name) - 1)) {
  for (j in (i + 1):length(group_name)) {
    group_ij = subset(sample_site, group %in% c(as.character(group_name[i]), as.character(group_name[j])))
    otu_ij = mymatrix[group_ij$names, ]
    adonis_result_otu_ij = adonis(otu_ij~group, group_ij, permutations = 999, distance = 'bray')
    res.temp = as.data.frame(adonis_result_otu_ij$aov.tab)[1,]
    rownames(res.temp) = paste(as.character(group_name[i]),'/',as.character(group_name[j]))
        result = rbind(result,res.temp)
  }
}
head(result,nrow(result))
