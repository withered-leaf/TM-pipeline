#eQTL analysis between host genotype and microbiome
#Written by Wei Guo
#Data:Nov 20,2020

#load library
library(MatrixEQTL)
## Location of the package with the data files.
#base.dir = find.package('MatrixEQTL');
# base.dir = '.';

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR ; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
#setwd("/storage/work/PancreaticCancerSubtypes/MatrixEQTL/out")
# Genotype file name
SNP_file_name = "immune_related_genotype_file";
# (Gene expression file name) Microbiome abundance file name
expression_file_name = "microbial_count_genus";
# Covariates file name
# Set to character() for no covariates
covariates_file_name = character();
# Output file name
output_file_name = tempfile();
# Only associations significant at this level will be saved
pvOutputThreshold = 1e-4;
# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

##Outliers in genotype. Minor allele frequency filtering.
maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]];
  maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)
cat('SNPs before filtering:',nrow(snps))
# snps$RowReorderSimple(maf>0.1);
snps$RowReorder(maf>0.1);
cat('SNPs after filtering:',nrow(snps))

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Outliers in expression. Quantile normalization
for( sl in 1:length(gene) ) {
  mat = gene[[sl]];
  mat = t(apply(mat, 1, rank, ties.method = "average"));
  mat = qnorm(mat / (ncol(gene)+1));
  gene[[sl]] = mat;
}
rm(sl, mat);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis

me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name);
## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
write.table(me$all$eqtls,"Results_LINEAR_Immune_microbiome",row.names=FALSE,sep="\t")
pdf(file="myplot_LINEAR_Immune_microbiome.qqplot.pdf")
plot(me,pch=16,cex=0.7)
dev.off()

#show(me$all$eqtls)