#microbial gene functional characterization
#Written by Wei Guo
#Data:Nov 20,2020

#load library
library(readxl)
library(matrixStats)
library(reshape)
library(ggplot2)
#
setwd("/data1/PDAC_files/code_ocean/")
##microbial genes were analysed based on eggNOG and GhostKOALA website, then processed by Perl script.
#load processed output.
KEGG_function<-read_excel("microbial_gene_functional_characterization.xlsx"
                          ,sheet = "Sheet1")
EGGNOG_function<-read_excel("microbial_gene_functional_characterization.xlsx",sheet = "Sheet2")
my_table<-KEGG_function
my_table<-my_table[,-1]
row.names(my_table)<-KEGG_function$KEGG
#calculate proportion
my_table<-as.matrix(my_table)
#all_geneset_number=229282
#basal_enriched_query_number=31905
sum_matrix<-matrix(rep(c(229282,31905),nrow(my_table)),nrow=nrow(my_table),byrow=T)
rela.table<-my_table/sum_matrix
#fisher exact test
size1 = sum_matrix[1,1]
size2 = sum_matrix[1,2]
fisher_exact<-function(x){
  row=c(x[1],x[2],size1-x[1],size2-x[2])
  alle<-matrix(row,nrow=2)
  p_value<-fisher.test(alle,alternative="two.sided")$p.value
  p_value}
pvalue<-apply(my_table,1,fisher_exact)
result<-cbind(my_table,pvalue)
final<-subset(result,pvalue<0.05)
#
signif<-ifelse(result[,3]<0.05,"*","")
the_data<-as.data.frame(rela.table)
colnames(the_data)[1]<-"Background"
the_data$signif<-signif
the_data$term<-row.names(the_data)
the_data<-melt(the_data,id.vars = c("term","signif"))
the_index<-c()
the_term=the_data$term[!duplicated(the_data$term)]
for(everyterm in the_term){
  myindex<-which(the_data$term==everyterm)
  myvalue=the_data$value[myindex]
  my_minvalue=min(myvalue)
  min_index=which(the_data$term==everyterm & the_data$value==my_minvalue)
  the_index<-c(the_index,min_index)
}

the_data$signif[the_index]<-""
res<-ggplot(the_data,aes(term,value,fill=variable))+geom_bar(stat="identity",position="dodge")+
  scale_fill_manual(values=c("cyan3","lightcoral"))+xlab("")+ylab("Proportion")+
  geom_text(aes(y=value+0.001,label=signif))+
  theme(legend.title = element_blank(),legend.position = c(0.9,0.9),
        axis.text.x = element_text(angle = 270,hjust = 0))
res
