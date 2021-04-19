#subtyping PDAC tumors using particular signatures
#Written by Wei Guo
#Data:Nov 20,2020

#load library
library(readxl)
library(ConsensusClusterPlus)
library(pheatmap)
#load RNA expression matrix 
setwd("/data")
rna.exp<-read.delim("PDA_rna_expression_matrix.csv",row.names = 1,stringsAsFactors = F)
#select gene signatures from Michelle publication
#pick top 25 gene signatures from each of sig2 sig10 sig1 sig6
gene.signatures<-read_excel("41588_2019_566_MOESM3_ESM.xlsx",sheet = "Supplementary Table 4",col_names = FALSE)
colnames(gene.signatures)<-gene.signatures[2,]
gene.signatures<-gene.signatures[c(-1,-2),]
Sig2<-gene.signatures$`Sig. 2 genes`[1:25]
Sig10<-gene.signatures$`Sig. 10 genes`[1:25]
Sig1<-gene.signatures$`Sig. 1 genes`[1:25]
Sig6<-gene.signatures$`Sig. 6 genes`[1:25]
#build the matrix
gene.exp<-rna.exp[c(Sig2,Sig10,Sig1,Sig6),]
gene.exp <- gene.exp[apply(gene.exp, 1, function(x){sum(is.na(x)) < ncol(gene.exp)/2}),]    
gene.exp<-as.matrix(gene.exp)
#build the row annotation
SigGroup<-c()
SigGroup[which(row.names(gene.exp) %in% Sig2)]<-"Sig2"
SigGroup[which(row.names(gene.exp) %in% Sig10)]<-"Sig10"
SigGroup[which(row.names(gene.exp) %in% Sig1)]<-"Sig1"
SigGroup[which(row.names(gene.exp) %in% Sig6)]<-"Sig6"
annotation_row<-data.frame(Signature=SigGroup,stringsAsFactors = T)
row.names(annotation_row)<-row.names(gene.exp)
##normalized by median
d = sweep(gene.exp,1, apply(gene.exp,1,median,na.rm=T))
##pheatmap 
##method: Average Hierarchical Clustering
##distance:Pearson Correlation
the_plot<-pheatmap(d,scale="row",color=colorRampPalette(c("blue","white","red"))(50),clustering_method="average",
                breaks = c(seq(-2,2,length=50)),show_rownames = F,show_colnames = F,main="",border_color = "NA",
                cluster_rows = F,clustering_distance_cols = "correlation",annotation_row = annotation_row,
                cutree_cols = 3,fontsize_row = 6,annotation_names_row = F,annotation_names_col = F)
colnames(d)[the_plot$tree_col$order]
##
my_group<-data.frame(Samples=colnames(d)[the_plot$tree_col$order],Subtypes=c(rep("Basal",17),rep("Hybrid",23),rep("Classical",22)))

  
save_pheatmap_pdf <- function(x, filename, width=10, height=6) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(the_plot, "heatmap.pdf")
###ConsensusCluster using the same clusterAlg and distance as pheatmap
res <- ConsensusClusterPlus(d, maxK = 5, reps = 1000, pItem = 0.8, pFeature = 1, clusterAlg = "hc", 
                            tmyPal=colorRampPalette(c("blue","white","red"))(50),distance="pearson",
                            title="Consensus_hc_4sig",innerLinkage="average", finalLinkage="average",
                            corUse = "complete.obs", seed=123456, plot='pdf', writeTable=T)
icl<-calcICL(res,title="ICL1",plot="pdf",writeTable=T)

###validate using signatures from Moffitt publication
gene.signatures_2<-read.delim("/data1/PDAC_files/RNA/Subtypes/Consensus_hc_T/25gene_for_Basal_Classical",header = F,row.names = NULL,stringsAsFactors = F)
the_signature2<-gene.signatures_2$V1
gene.exp<-rna.exp[the_signature2,]
gene.exp <- gene.exp[apply(gene.exp, 1, function(x){sum(is.na(x)) < ncol(gene.exp)/2}),]
##exclude an abnormal signature in our cohort 
gene.exp<-gene.exp[!row.names(gene.exp) %in% c("LEMD1"),]
gene.exp<-as.matrix(gene.exp)
#build the row annotation
SigGroup<-c()
SigGroup[which(row.names(gene.exp) %in% the_signature2[1:25])]<-"Basal_like"
SigGroup[which(row.names(gene.exp) %in% the_signature2[26:50])]<-"Classical"
annotation_row<-data.frame(Signature=SigGroup,stringsAsFactors = T)
row.names(annotation_row)<-row.names(gene.exp)

##normalized by median
d = sweep(gene.exp,1, apply(gene.exp,1,median,na.rm=T))
##pheatmap 
##method: Average Hierarchical Clustering
##distance:Pearson Correlation
the_plot<-pheatmap(d,scale="row",color=colorRampPalette(c("blue","white","red"))(50),clustering_method="average",
                   breaks = c(seq(-2,2,length=50)),show_rownames = T,show_colnames = F,main="",border_color = "NA",
                   cluster_rows = T,clustering_distance_cols = "correlation",annotation_row = annotation_row,
                   cutree_cols = 3,fontsize_row = 6,annotation_names_row = F,annotation_names_col = F)
colnames(d)[the_plot$tree_col$order]
my_group2<-data.frame(Samples=colnames(d)[the_plot$tree_col$order],Subtypes=c(rep("Classical",27),rep("Basal",17),rep("Others",18)))
#
save_pheatmap_pdf(the_plot, "heatmap2.pdf")
###ConsensusCluster using the same clusterAlg and distance as pheatmap
res <- ConsensusClusterPlus(d, maxK = 5, reps = 1000, pItem = 0.8, pFeature = 1, clusterAlg = "hc", 
                            tmyPal=colorRampPalette(c("blue","white","red"))(50),distance="pearson",
                            title="Consensus_hc_2sig",innerLinkage="average", finalLinkage="average",
                            corUse = "complete.obs", seed=123456, plot='pdf', writeTable=T)
icl<-calcICL(res,title="ICL2",plot="pdf",writeTable=T)
##
###sample alignments
dd<-merge(my_group,my_group2,by="Samples",all=T,sort=F)
colnames(dd)[2:3]<-c("Michelle et al.","Moffitt et al.")
row.names(dd)<-dd$Samples
dd<-dd[,-1]
dd<-dd[,c(2,1)]
annotation_col<-dd
anno_colors=list(`Michelle et al.`=c(Classical="lightskyblue",Basal="lightcoral",Hybrid="azure4"),
                 `Moffitt et al.`=c(Classical="lightskyblue",Basal="lightcoral",Others="lightgray"))
my_matrix<-rna.exp[1,]
my_matrix<-my_matrix[,row.names(dd)]
the_plot2<-pheatmap(my_matrix,scale="row",color=colorRampPalette(c("white","white","white"))(50),clustering_method="average",
                    breaks = c(seq(-2,2,length=50)),show_rownames = F,show_colnames = F,main="",border_color = "white",
                    cluster_rows = F,cluster_cols = F,clustering_distance_cols = "correlation",annotation_row = NULL,annotation_col = annotation_col,
                    fontsize_row = 6,annotation_names_row = F,annotation_names_col = F,cellwidth = 6,cellheight = 12,legend = F,annotation_colors = anno_colors)
