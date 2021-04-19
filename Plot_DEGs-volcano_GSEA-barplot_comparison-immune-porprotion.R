#presentation of differential gene expression, gene set enrichment analysis and Cibersortx.
#Written by Wei Guo
#Data:Nov 20,2020

setwd("/data")
##volcano plot
library(ggplot2)
library(ggrepel)
library(readxl)
library(reshape)
library(ggpubr)
gene.table<-read.delim("BasalvsClassical_DEGs_all",header = T,row.names = 1,stringsAsFactors = F)
data=as.data.frame(gene.table)
data$threshold <- as.factor(ifelse(data$adj.P.Val < 0.01 & abs(data$logFC) > 1,ifelse(data$logFC > 1 ,'Up','Down'),'Not'))
data$threshold[is.na(data$threshold)]=c('Not')
p=ggplot(data=data,aes(x=logFC, y =-log10(adj.P.Val),colour=threshold,fill=threshold)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_point(alpha=0.8, size=1.2)+
  xlim(c(-6, 6)) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6)+
  geom_hline(yintercept = -log10(0.01),lty=4,col="grey",lwd=0.6)+
  theme(legend.position="right",
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "Times", size=8),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=12),
        axis.text.y = element_text(face="bold",  color="black", size=12),
        axis.title.x = element_text(face="bold", color="black", size=12),
        axis.title.y = element_text(face="bold",color="black", size=12))+
  labs(x="log2 (fold change)",y="-log10 (p-value)",title="")
p
pdf('volcano.pdf',width = 8,height = 6)

#
####barplot
my_nes<-read_excel("NES.xlsx",sheet = "Sheet1")
my_nes$ID<-sub("HALLMARK_","",my_nes$ID)
my_nes<-my_nes[order(my_nes$NES,decreasing = T),]
res<-ggplot(data = my_nes,
            aes(x = reorder(ID,NES), y = NES,
                fill = NES > 0))+  geom_bar(stat = "identity")+  coord_flip()+guides(fill=F)+
  xlab("")+ylab("Normalized enrichment score")+scale_fill_manual(values=c("deepskyblue","lightcoral"))+
  geom_text(aes(y=0,label = ID,hjust = ifelse(NES >= 0, 1.05, -0.05)),size=3,color="gray20")+
  theme_bw()+theme(axis.text.y=element_blank(),axis.ticks.y =element_blank(),plot.background = element_blank(),
                   panel.border = element_blank(),panel.background = element_blank())

res
pdf(file="HALLMARK.pdf",width = 8,height = 6)

##CibersortX
CIBER_result<-read.csv("CIBERSORTx_Job6_Results.csv",row.names=1,stringsAsFactors=FALSE)
mydata<-data.frame(sample=row.names(CIBER_result),CIBER_result)
del=which(colnames(mydata) %in% c("P.value","Correlation","RMSE"))
mydata<-mydata[,-del]
#Groupinfo
Groupinfo<-read.delim("Sample_group.txt",stringsAsFactors = F)
Groupinfo$Subtype<-factor(Groupinfo$Subtype,levels = c("Classical","Hybrid","Basal"))
##
mydata<-mydata[Groupinfo$Samples,]
mydata$group<-Groupinfo$Subtype
mydata<-melt(mydata,id.vars = c("sample","group"))
##subset
mydata<-subset(mydata,variable %in% c("B.cells.memory","T.cells.follicular.helper","NK.cells.resting","Macrophages.M1",
                                      "Mast.cells.activated"))
kruskal_res<-compare_means(value ~ group, data = mydata, group.by = "variable",method = "kruskal.test")
res1<-ggplot(mydata,aes(variable,value,fill=group))+geom_boxplot(aes(fill=group),outlier.color=NA)+
    scale_fill_manual(values = c("deepskyblue","gray40","brown1"))+xlab("")+ylab("Proportion")+coord_flip()+
  theme(legend.title = element_blank(),legend.position = "bottom")+
  stat_compare_means(aes(group = group), label = "p.signif",hide.ns=T)
res1
ggsave("immune_boxplot.pdf",res1,width=5,height = 3)
