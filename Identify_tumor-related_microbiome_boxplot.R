#Identify differential abundant microorganisms 
#Written by Wei Guo
#Data:Nov 20,2020

#load library
library(ggplot2)
library(ggrepel)
library(dplyr)
library(pheatmap)
library(ggpubr)
#load microbial abundance matrix from kraken2+Bracken output.
setwd("/data1/PDAC_files/code_ocean/")
#the_taxonomic_level="genus"
myfile="PDAC_samples_count_genus"
bac.exp<-read.delim(file = myfile,row.names = 1,stringsAsFactors = F)
#
#### load group information
Groupinfo<-read.delim("Sample_group.txt",stringsAsFactors = F)
Groupinfo$Subtype<-factor(Groupinfo$Subtype,levels = c("Classical","Hybrid","Basal"))
bac.exp<-bac.exp[,Groupinfo$Samples]
#
group1=which(Groupinfo$Subtype=="Classical")
group2=which(Groupinfo$Subtype=="Hybrid")
group3=which(Groupinfo$Subtype=="Basal")
wilcox_test<-function(x){
  p1_value<-wilcox.test(x[group1],x[group2],paired=F)$p.value
  p2_value<-wilcox.test(x[group2],x[group3],paired=F)$p.value
  p3_value<-wilcox.test(x[group1],x[group3],paired=F)$p.value
  p_out<-c(p1_value,p2_value,p3_value)
  p_out}
wil_pvalue<-apply(bac.exp,1,wilcox_test)
wil_pvalue<-t(wil_pvalue)
wil_pvalue<-round(wil_pvalue,digit=5)
colnames(wil_pvalue)<-c("Classical.vs.Hybrid","Hybrid.vs.Basal","Classical.vs.Basal")

kru_test<-function(x){
  my_data<-data.frame(value=x,groups=Groupinfo$Subtype)
  p<-kruskal.test(value~groups,my_data)$p.value
  return(p)
}
kru_pvalue<-apply(bac.exp, 1, kru_test)

tend<-function(x){
  fc1<-log2(round(mean(x[group2])/mean(x[group1]),digit=2))
  fc2<-log2(round(mean(x[group3])/mean(x[group2]),digit=2))
  fc3<-log2(round(mean(x[group3])/mean(x[group1]),digit=2))
  fc_out<-c(fc1,fc2,fc3)
  return(fc_out)}

FC<-apply(bac.exp,1,tend)
FC<-t(FC)
colnames(FC)<-c("FC_Classical.vs.Hybrid","FC_Hybrid.vs.Basal","FC_Classical.vs.Basal")
fdr_pvalue<-p.adjust(kru_pvalue,method="fdr")
data<-cbind(FC,wil_pvalue,kru_pvalue,fdr_pvalue)
# filter which pvalue <0.05
final<-as.data.frame(data)[which(kru_pvalue<0.05),]
final<-final[order(final$kru_pvalue),]

#volcano plot
data<-as.data.frame(data)
data$threshold<-as.factor(ifelse(data$Classical.vs.Basal < 0.05   ,ifelse(data$FC_Classical.vs.Basal > 0, "Enriched","Depleted"),"Not"))
data$microname=row.names(data)
p=ggplot(data=data,aes(x=FC_Classical.vs.Basal, y =-log10(Classical.vs.Basal),colour=threshold,fill=threshold)) +
  scale_color_manual(values=c( "blue","red","gray"))+
  geom_point(alpha=0.8, size=1.2)+
  xlim(c(-6, 6)) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept=0,lty=4,col="grey",lwd=0.6)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6)+
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
#p+geom_text_repel(data=filter(data, Classical.vs.Basal<0.05&FC_Classical.vs.Basal>0), aes(label=microname),show.legend=F,cex=2.5)

##heatmap
sig.bac.list<-row.names(data[which(data$kru_pvalue<0.05),])
my.sig.bac.exp<-bac.exp[sig.bac.list,]
annotation_col<-data.frame(Subtypes=Groupinfo$Subtype,stringsAsFactors = F)
row.names(annotation_col)=colnames(my.sig.bac.exp)
anno_colors=list(Subtypes=c(Classical="deepskyblue",Hybrid="gray20",Basal="brown1"))
the_plot<-pheatmap(my.sig.bac.exp,scale="row",annotation_col = annotation_col,breaks = c(seq(-2,2,length=100)),show_rownames = T,
                   show_colnames = F,main="",border_color = "white",annotation_names_col = F,fontsize_row = 4,
                   annotation_colors = anno_colors,cluster_cols = T)

save_pheatmap_pdf <- function(x, filename, width=10, height=8) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(the_plot, "heatmap.pdf")

#sig.boxplot for each candidate microbe
candidate.microbes<-c("Acinetobacter","Pseudomonas","Sphingopyxis")
for(microname in candidate.microbes)
{
  test_cycle<-as.data.frame(t(my.sig.bac.exp[row.names(my.sig.bac.exp)==microname,]))
  names(test_cycle)="counts"
  test_cycle$groups=Groupinfo$Subtype
  filename<-paste(microname,"counts",sep="_")
  #svg(file=filename)
  ggplot(test_cycle,aes(groups,counts))+geom_boxplot(aes(color=groups),outlier.color = NA)+geom_jitter(aes(color=groups),size=0.8)+
    scale_color_manual(values=c("deepskyblue","gray20","brown1"))+ggtitle(microname)+theme(plot.title = element_text(hjust = 0.5),legend.position = "none")+xlab("")+
    stat_compare_means()
  #dev.off()
  ggsave(paste(filename,"pdf",sep="."),device = pdf,width = 4.5,height = 4.5)
}
