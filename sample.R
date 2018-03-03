# Load required packages
library("Biobase"); # for eset
library("genefilter"); # for rowttests
library("gtools"); # for foldchange, foldchange2logratio
library("pheatmap"); # for pheatmap

base.means.and.fold.change <- function(x, feature_indices) {
	this.result <- matrix(NA, nrow = 1, ncol = 5);
	colnames(this.result) <- c("baseMean", "baseMeanA", "baseMeanB", "foldChange", "log2FoldChange");
	this.result[1] <- mean(x);
	this.result[2] <- mean(x[feature_indices[[1]]]);
	this.result[3] <- mean(x[feature_indices[[2]]]);
	this.result[4] <- foldchange(this.result[,3], this.result[,2]);
	this.result[5] <- foldchange2logratio(this.result[,4]);
	return(this.result)
}

# Normalized values here
conds <- array(c(), dim=c());

eset <- new("ExpressionSet", exprs=data.matrix(norm_rna));
ttest <- rowttests(exprs(eset), factor(conds));
# Map features to columns by indices
features <- sort(unique(conds));
feature_indices <- list();
for (i in (1:length(features))) {
	feature_indices[[i]] <- which(conds == features[i]  %in% c(TRUE))
}
# Calculate mean, mean by feature, fold change, and log_2 ratio of fold change
base.mean.and.fold.change <- apply(norm_rna, 1, base.means.and.fold.change, feature_indices = feature_indices);
base.mean.and.fold.change <- t(base.mean.and.fold.change);
colnames(base.mean.and.fold.change) <- c("baseMean", "baseMeanA", "baseMeanB", "foldChange", "log2FoldChange");
# Store calculations into ttest data frame
ttest <- cbind(ttest, base.mean.and.fold.change);
# Calculate adjusted p-value
ttest$p.value.adj <- p.adjust(ttest$p.value, method = "BH");
# Add id column, change rownames to numbers, reorder columns
ttest$id <- rownames(ttest);
rownames(ttest) <- 1:nrow(ttest);
ttest <- ttest[,c(10,4:8,3,9,1:2)];
sort.ttest <- ttest[order(ttest$p.value.adj),]; # sort results by p-values

pval_row_index <- which((ttest$p.value.adj < 0.05) %in% c(TRUE));
filtered_norm_rna <- norm_rna[pval_row_index, ]; # get normalized data with p.value.adj < 0.05

# Convert _ to whitespace
col_names <- colnames(data.matrix(filtered_norm_rna));
col_names <- gsub("(sample_number__)", " ", col_names);
col_names <- gsub("_", " ", col_names);
colnames(filtered_norm_rna) <- col_names;
# Adjust the margins
row_names <- rownames(data.matrix(filtered_norm_rna));
max_char_row <- max(nchar(row_names));
max_char_col <- max(nchar(col_names));
# Adjust image size
# size for plot + size for label margins + size for dendrogram in pixels + outer margins
heatmap_height <- length(row_names) * 30 + max_char_col * 8 + 50 + 50;
heatmap_width <- length(col_names) * 30 + max_char_row * 8 + 50 + 50 + 75; # add 75 for legend
png("heatmap.png", height = heatmap_height, width = heatmap_width, res = 72);
pheatmap(data.matrix(filtered_norm_rna), cluster_cols = TRUE, legend = FALSE, scale = "row", border_color = NA, cellwidth = 30, cellheight = 30, treeheight_row = 75, treeheight_col = 75, annotation_legend = FALSE, fontsize = 12, fontfamily = "mono", fontface = "plain", color = colorRampPalette(c("green", "black", "red"))(100));
dev.off();