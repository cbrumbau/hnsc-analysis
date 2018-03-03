library('TCGAbiolinks'); # For data retrieval and analysis
library('SummarizedExperiment') # For processing SummarizedExperiment objects
library('igraph') # Plot the Cox-regression analysis
#library('CombinePValue'); # For Fisher method of analysis (if multiple platform results for meta analysis)

# Query for RNASeqv2 data in HNSC
query <- TCGAquery(tumor = c('HNSC'), level = 3, platform = c('IlluminaHiSeq_RNASeqV2', 'IlluminaGA_RNASeqV2', 'IlluminaHiSeq_mRNA_DGE', 'IlluminaGA_mRNA_DGE'))
# Download the query results
TCGAdownload(query, path = './dataHNSC', type = 'rsem.genes.results')
# Read in the data
prep.df <- TCGAprepare(query, './dataHNSC', type = 'rsem.genes.results')
# Optionally look at the preprocessed data
#TCGAanalyze_Preprocessing(prep.df, file="preprocessing.png", width = 1000, height = 1000)
# Normalize
raw_counts.mat <- assay(prep.df, "raw_counts")
data_norm.mat <- TCGAanalyze_Normalization(tabDF = raw_counts.mat, geneInfo =  geneInfo)
# Filter background noise by quantile
data_filter.mat <- TCGAanalyze_Filtering(tabDF = data_norm.mat, method = "quantile", qnt.cut = 0.25)
# Get the conditions for samples
normal.list <- TCGAquery_SampleTypes(barcode = colnames(data_filter.mat), typesample = c('NT'))
tumor.list <- TCGAquery_SampleTypes(barcode = colnames(data_filter.mat), typesample = c('TP'))
# Use edgeR for differential expression analysis
data_diffexp.tab <- TCGAanalyze_DEA(mat1 = data_filter.mat[,normal.list], mat2 = data_filter.mat[,tumor.list], Cond1type = "Normal", Cond2type = "Tumor", fdr.cut = 0.01, logFC.cut = 1, method = "glmLRT")
# Generate table of differentially expressed mRNAs with expression levels
data_diffexp_filter.tab <- TCGAanalyze_LevelTab(data_diffexp.tab, "Tumor", "Normal", data_filter.mat[,tumor.list], data_filter.mat[,normal.list])
# Get the clinical data
data_clinical.df <- TCGAquery_clinic(tumor = c('HNSC'), clinical_data_type = "clinical_patient")
# Perform Kaplan-Meier analysis
data_survival.tab <- TCGAanalyze_SurvivalKM(clinical_patient = data_clinical.df, dataGE = data_filter.mat, Genelist = rownames(data_diffexp.tab), Survresult = FALSE, ThreshTop = 0.67, ThreshDown = 0.33, p.cut = 0.05)
# Perform Cox-regression analysis
download.file('http://dnet.r-forge.r-project.org/RData/1.0.7/org.Hs.string.RData', './org.Hs.string.RData', 'auto')
load('org.Hs.string.RData')
coxnet.igraph <- TCGAvisualize_SurvivalCoxNET(data_clinical.df, data_filter.mat, Genelist = rownames(data_survival.tab), scoreConfidence = 700, org.Hs.string = org.Hs.string, titlePlot = "Cox regression analysis of HNSC")
# Graph the results
plot.igraph(coxnet.igraph, vertex.size=3,vertex.label=NA, layout=layout.fruchterman.reingold(g, niter=10000, area=30*vcount(g)^2))
# Get the differentially expressed genes
gene.list <- rownames(data_diffexp_filter.tab)
# Do enrichment analysis of genes
result_enrich.tab <- TCGAanalyze_EAcomplete(TFname="HNSC Normal vs Tumor", gene.list)
# Plot the enrichment analysis results
TCGAvisualize_EAbarplot(tf = rownames(result_enrich.tab$ResBP), GOBPTab = result_enrich.tab$ResBP, GOCCTab = result_enrich.tab$ResCC, GOMFTab = result_enrich.tab$ResMF, PathTab = result_enrich.tab$ResPat, nRGTab = gene.list, nBar = 20, filename='HNSC_GO_Analysis.pdf')
# Generate volcano plot of differential expression analysis
TCGAVisualize_volcano(data_diffexp.tab$logFC, data_diffexp.tab$FDR, filename = "HNSC_volcano.png", x.cut = 5, y.cut = 10^-50, names = rownames(data_diffexp.tab), names.size = 2, xlab = " Gene expression fold change (Log2)", legend = "State", title = "Volcano plot (Normal vs Tumor)", width = 20)
# Plot a heatmap of the DEGs
TCGAvisualize_Heatmap(t(data_filter.mat), col.metadata =  clin_subt[,c("bcr_patient_barcode", "groupsHC", "histological_type", "IDH.codel.subtype")], col.colors =  list(groupsHC = c("EC1"="black", "EC2"="red", "EC3"="blue", "EC4"="green3"), histological_type=c("Astrocytoma"="navy", "Oligoastrocytoma"="green3", "Oligodendroglioma"="red"), IDH.codel.subtype = c("IDHmut-codel"="tomato", "IDHmut-non-codel"="navy", "IDHwt"="gold","NA"="white")), sortCol = "groupsHC", type = "expression", scale = "row", title = "Heatmap from concensus cluster", cluster_rows = TRUE)
# Use PCA to visualize clusters
#TCGAvisualize_PCA(data_filter.mat, data_diffexp_filter.tab, ntopgenes = 200)