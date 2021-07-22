#Correlation between host gene expression and microbiome  
#Written by Wei Guo
#Data:Nov 20,2020

#load library
library(clusterProfiler,quietly = T)
library(MEGENA,quietly = T)
library(psych,quietly = T)
library(reshape2,quietly = T)
library(readxl)
library(msigdf)
library(dplyr)
library(org.Hs.eg.db)
library(GSEABase)
library(GSVA)

##load input data. microbiome 
bac.exp<-read.delim(file = "PDAC_samples_count_genus",row.names = 1,stringsAsFactors = F)
zero.count=apply(bac.exp,1,function(x){length(which(x>0))/length(x)})
bac.exp=bac.exp[which(zero.count>0.2),]
#load human gene expression
rna.exp<-read.delim("PDA_rna_expression_matrix.csv",row.names = 1,stringsAsFactors = F)
gene.table<-read.delim("BasalvsClassical_DEGs_V2",header = T,row.names = 1,stringsAsFactors = F)
deg.gene.name<-row.names(gene.table)
human.deg.exp=rna.exp[deg.gene.name,]
#
#kegg pathway enrichment
eg = bitr(deg.gene.name, fromType="ENSEMBL", toType="ENTREZID",OrgDb = "org.Hs.eg.db")
kk <- enrichKEGG(gene         = eg$ENTREZID,
                 organism     = 'hsa',
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.5
)
kk=setReadable(kk, OrgDb = org.Hs.eg.db,keyType ="ENTREZID")
View(as.data.frame(kk))
#data.frame(kk)
pathway=data.frame(kk)
#GO term enrichment
##
ego<-enrichGO(OrgDb="org.Hs.eg.db", 
              gene = deg.gene.name,
              keyType = "SYMBOL",
              ont = "ALL",
              pvalueCutoff = 0.05,
              pAdjustMethod = "BH",
              qvalueCutoff = 0.2
              
)
View(as.data.frame(ego))
#
orderedgene<-sigres[order(sigres$logFC,decreasing = T),]
GeneList<-orderedgene[,1]
names(GeneList)<-row.names(orderedgene)
##hallmark enrichment analysis
thegeneset<-msigdf.human %>% 
  filter(category_code == "hallmark") %>% 
  select(geneset, symbol) %>% as.data.frame

x <- enricher(names(GeneList), TERM2GENE = thegeneset,pvalueCutoff = 0.05,pAdjustMethod = "none",qvalueCutoff = 0.5)
View(as.data.frame(x))
##
#PCA -- get PC1
get.PC1=function(data){
  pca=prcomp(data,scale= T)
  pc1=pca$rotation[,1]
}
psych_cor <- function(mat.1,mat.2,SigThres=0.05){
  corres <- corr.test(t(mat.1),t(mat.2),ci=F,use='pairwise',method='pearson',adjust=
                        'none')
  cr <- corres$r
  cp <- corres$p
  cres <- corr.test(t(mat.1),t(mat.2),ci=F,use='pairwise',method='pearson',adjust='BH')
  cpadj <- cres$p
  mcor <- melt(cr)
  mpval <- melt(cp)
  mpadj <- melt(cpadj)
  tb <- cbind(mcor[which(mpval$value<=SigThres),],mpval$value[which(mpval$value<=SigThres)])
  colnames(tb) <- c('feature1','feature2','correlation','p_value')
  return(tb)
}
#
pathways.PC1s=sapply(1:dim(selected_pathways)[1],function(n){
  gene.string=selected_pathways$geneID[n]
  genes=unlist(strsplit(gene.string,split="/"))
  my.exp=human.exp[row.names(human.exp) %in% genes,]
  pc1=get.PC1(my.exp)
  return(pc1)
})
colnames(pathways.PC1s)=paste(selected_pathways$ID,selected_pathways$Description)
#
GOs.PC1s=sapply(1:dim(selected_GOs)[1],function(n){
  gene.string=selected_GOs$geneID[n]
  genes=unlist(strsplit(gene.string,split="/"))
  my.exp=human.exp[row.names(human.exp) %in% genes,]
  pc1=get.PC1(my.exp)
  return(pc1)
})
colnames(GOs.PC1s)=paste(selected_GOs$ID,selected_GOs$Description)
##
hallmark.PC1s=sapply(1:dim(selected_hallmark)[1],function(n){
gene.string=selected_hallmark$geneID[n]
genes=unlist(strsplit(gene.string,split="/"))
my.exp=human.exp[row.names(human.exp) %in% genes,]
pc1=get.PC1(my.exp)
return(pc1)
})
colnames(hallmark.PC1s)=selected_hallmark$Description
#
hallset<-getGmt("h.all.v7.1.symbols.gmt")
hallES<-gsva(as.matrix(human.exp),gset.idx.list = hallset,kcdf="Poisson",parallel.sz=4)#or Gaussian
row.names(hallES)<-sub("HALLMARK_","",row.names(hallES))
keggset<-getGmt("c2.cp.kegg.v7.1.symbols.gmt")
keggES<-gsva(as.matrix(human.exp),gset.idx.list = keggset,kcdf="Poisson",parallel.sz=4)
row.names(keggES)<-sub("KEGG_","",row.names(keggES))
#
cor.bac.pathway.res=psych_cor(bac.exp,t(pathways.PC1s))
cor.bac.pathway.all=psych_cor(bac.exp,t(pathways.PC1s),SigThres=1)
#
cor.bac.GO.res=psych_cor(bac.exp,t(GOs.PC1s))
cor.bac.GO.all=psych_cor(bac.exp,t(GOs.PC1s),SigThres = 1)
#
cor.bac.hallmark.res=psych_cor(bac.exp,t(hallmark.PC1s))
cor.bac.hallmark.all=psych_cor(bac.exp,t(hallmark.PC1s),SigThres = 1)
#
cor.bac.gsva.hallmark.res=psych_cor(bac.exp,hallES)
cor.bac.gsva.hallmark.all=psych_cor(bac.exp,hallES,SigThres = 1)
#
cor.bac.gsva.kegg.res=psych_cor(bac.exp,keggES)
cor.bac.gsva.kegg.all=psych_cor(bac.exp,keggES,SigThres = 1)
#basal-like tumor enriched bacteria
enriched.sig.bac<-c("Massilia","Corynebacterium","Brevundimonas","Paracoccus","Sphingomonas","Mitsuaria","Streptomyces",
                    "Brachybacterium","Polaromonas","Mesorhizobium","Blastomonas","Pseudomonas","Caulobacter",
                    "Sphingopyxis","Collimonas","Ensifer","Laribacter","Tessaracoccus","Shinella","Acinetobacter",
                    "Bosea","Chelatococcus","Comamonas","Flavobacterium","Methyloversatilis","Melaminivora",
                    "Ochrobactrum","Rhodobacter","Agrobacterium","Thiomonas","Paucibacter",
                    "Rhodopseudomonas","Sinorhizobium","Bordetella","Altererythrobacter","Alicycliphilus",
                    "Sphingorhabdus","Porphyrobacter","Methylibium","Rubrivivax","Novosphingobium",
                    "Erythrobacter","Nocardioides")
### merge process,  
#Gsva c2.kegg.pathway all.ES(human.exp) ~ bac.exp
my_pickup<-c("PANCREATIC_CANCER","DNA_REPLICATION")
c2.kegg.gsva.all<-cor.bac.gsva.kegg.all[which(cor.bac.gsva.kegg.all$feature2 %in% my_pickup),]
c2.kegg.gsva.res<-cor.bac.gsva.kegg.res[which(cor.bac.gsva.kegg.res$feature1 %in% enriched.sig.bac),]
c2.kegg.gsva.res<-c2.kegg.gsva.res[which(c2.kegg.gsva.res$feature2 %in% my_pickup),]
#hallmark  Gsva hallmark.all  all.ES(human.exp) ~ bac.exp
my_pickup<-c("WNT_BETA_CATENIN_SIGNALING","G2M_CHECKPOINT","PI3K_AKT_MTOR_SIGNALING","E2F_TARGETS"
             ,"P53_PATHWAY","BILE_ACID_METABOLISM","PANCREAS_BETA_CELLS")
hallmark.gsva.all<-cor.bac.gsva.hallmark.all[which(cor.bac.gsva.hallmark.all$feature2 %in% my_pickup),]
hallmark.gsva.res<-cor.bac.gsva.hallmark.res[which(cor.bac.gsva.hallmark.res$feature1 %in% enriched.sig.bac),]
hallmark.gsva.res<-hallmark.gsva.res[which(hallmark.gsva.res$feature2 %in% my_pickup),]
## selected.KEGG.pathways
my_pickup<-c("hsa04911 Insulin secretion","hsa05204 Chemical carcinogenesis","hsa04972 Pancreatic secretion",
             "hsa04015 Rap1 signaling pathway","hsa04014 Ras signaling pathway","hsa04010 MAPK signaling pathway")
selected.kegg.cor.all<-cor.bac.pathway.all
selected.kegg.cor.res<-cor.bac.pathway.res[which(cor.bac.pathway.res$feature1 %in% enriched.sig.bac),]
selected.kegg.cor.res<-selected.kegg.cor.res[which(selected.kegg.cor.res$feature2 %in% my_pickup),]
### selected.hallmark categories
my_pickup<-c("KRAS_SIGNALING_UP","EPITHELIAL_MESENCHYMAL_TRANSITION","COMPLEMENT")
selected.hall.cor.all<-cor.bac.hallmark.all[which(cor.bac.hallmark.all$feature2 %in% my_pickup),]
selected.hall.cor.res<-cor.bac.hallmark.res[which(cor.bac.hallmark.res$feature1 %in% enriched.sig.bac),]
selected.hall.cor.res<-selected.hall.cor.res[which(selected.hall.cor.res$feature2 %in% my_pickup),]
###selected.GOs terms
my_pickup<-c("GO:0009913 epidermal cell differentiation","GO:0002237 response to molecule of bacterial origin",
             "GO:0032496 response to lipopolysaccharide")
selected.GOs.cor.all<-cor.bac.GO.all[which(cor.bac.GO.all$features %in% my_pickup)]
selected.GOs.cor.res<-cor.bac.GO.res[which(cor.bac.GO.res$feature1 %in% enriched.sig.bac),]
selected.GOs.cor.res<-selected.GOs.cor.res[which(selected.GOs.cor.res$feature2 %in% my_pickup),]
####
final.selected.cor.res<-rbind(c2.kegg.gsva.res,hallmark.gsva.res,selected.kegg.cor.res,selected.hall.cor.res,selected.GOs.cor.res)
final.selected.cor.all<-rbind(c2.kegg.gsva.all,hallmark.gsva.all,selected.kegg.cor.all,selected.hall.cor.all,selected.GOs.cor.all)
## genernal plot parameter
calculate.r<-function(cor.matrix.all,cor.matrix.res){
  cor.all.r<-cor.matrix.all[,1:3]
  cor.all.r.mtx<-dcast(cor.all.r,feature1~feature2)
  cor.input.res<-cor.matrix.res[which(cor.matrix.res$feature1 %in% enriched.sig.bac),]
  my.uniq.bac.list<-as.character(cor.input.res$feature1[!duplicated(cor.input.res$feature1)])
  my.uniq.term.list<-as.character(cor.input.res$feature2[!duplicated(cor.input.res$feature2)])
  cor.r.candidates<-cor.all.r.mtx[which(cor.all.r.mtx$feature1 %in% my.uniq.bac.list),c(1,which(colnames(cor.all.r.mtx) %in% my.uniq.term.list))]
  row.names(cor.r.candidates)<-cor.r.candidates$feature1
  cor.r.candidates<-cor.r.candidates[,-1]
  return(cor.r.candidates)
}
calculate.p<-function(cor.matrix.all,cor.matrix.res){
  cor.all.p<-cor.matrix.all[,c(1,2,4)]
  cor.all.p.mtx<-dcast(cor.all.p,feature1~feature2)
  cor.input.res<-cor.matrix.res[which(cor.matrix.res$feature1 %in% enriched.sig.bac),]
  my.uniq.bac.list<-as.character(cor.input.res$feature1[!duplicated(cor.input.res$feature1)])
  my.uniq.term.list<-as.character(cor.input.res$feature2[!duplicated(cor.input.res$feature2)])
  cor.p.candidates<-cor.all.p.mtx[which(cor.all.p.mtx$feature1 %in% my.uniq.bac.list),c(1,which(colnames(cor.all.p.mtx) %in% my.uniq.term.list))]
  row.names(cor.p.candidates)<-cor.p.candidates$feature1
  cor.p.candidates<-cor.p.candidates[,-1]
  return(cor.p.candidates)
}
the.cor.r.candidates<-calculate.r(final.selected.cor.all,final.selected.cor.res)
the.cor.p.candidates<-calculate.p(final.selected.cor.all,final.selected.cor.res)
pheatmap(the.cor.r.candidates,color = colorRampPalette(c("lightskyblue","white","hotpink"))(50),
         border_color = "white",display_numbers = matrix(ifelse(the.cor.p.candidates <= 0.05, ifelse(the.cor.p.candidates <=0.01, "*","+"), ""), nrow(the.cor.p.candidates)),
         fontsize_number = 5,fontsize=6.5,fontsize_row = 6.5)

