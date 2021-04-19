# Tumor microbiome analysis pipeline
R codes to analyze the microbial communities from pancreatic tumors sequencing. The pipeline includes tumors clustering, microbial profiles building, differential abundant microorganisms identification, microbial diversity measurement, microbial functional characterization, correlation between host gene expression and microbiome, association between host genetics and microbiome etc. 

The code files here are linked to the work "Tumor microbiome contributes to an aggressive phenotype in the basal-like subtype of pancreatic cancer" by Wei Guo et al. 

# Contents

* R codes to perform the ConsensusCluster and plot the heatmap of gene expression based on given signatures.

`PDAC_tumors_clustering_plot_heatmap.R`

* R codes to present the volcano plot of differentially expressed genes between tumor subtypes, to plot the barplot of Gene Set Enrichment Analysis (GSEA) result, and to show the difference of immune cell infiltration between tumor subtypes.

`Plot_DEGs-volcano_GSEA-barplot_comparison-immune-porprotion.R`

* R codes to depict the predominant microbiota composition landscape of PDAC tumors. Here we take an example: in taxonomic level of genus. The codes involve the measurement of microbial diversity using richness (the number of observed taxonomic units), Shannon index of alpha-diversity, and Bray-Curtis metric distance of beta-diversity. In addition, principal coordinates analysis (PCoA) was carried out to show the difference between tumor subtypes.

`Microbial_profile_diversity_measure_PCoA_plot.R`

* R codes to identify the significantly differential abundant microorganisms between groups using Kruskal-Wallis test. Plot Boxplots to show the abundance level of candidate bacteria between tumor subtypes.

`Identify_tumor-related_microbiome_boxplot.R`

* R codes to compare the microbial gene functional characterization results from eggNOG and GhostKOALA website, and to plot the barplot of the OG functional categories and level 2 of KEGG functional categories.

`Microbial_genes_functional_characteriazation.R`

* R codes to carry out the correlation analysis between host gene expression and microbiome. In brief, we performed gene enrichment analysis against several functional categories, including KEGG pathways, GO terms and the hallmark gene set, using a hypergeometric test by clusterProfiler. The first principal component (PC1) was estimate to represent the general expression level of the functional module. The correlation between functional module and microbial abundance was analyzed using the Pearson algorithm in the R package psych. 

`Correlation_between_host_gene_expression_and_microbiome.R`

* R codes to perform the association analysis between host genetics and microbiome. To link the microbial genera to genetic variation, we treated the abundance of genera as quantitative traits, and quantitative trait locus (QTL) mapping was carried out using the R package MatrixEQTL. We chose a linear model based on the assumption that genotypes have only additive effects on microbial abundance.

`eQTL_analysis_between_host_genotype_and_microbiome.R`

The example files can be found in the data folder. 
