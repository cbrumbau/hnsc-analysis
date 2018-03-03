# Load required packages
library("Biobase"); # for eset
library("GEOquery"); # for getGEO
library("pheatmap"); # for pheatmap

anova.matrix <- function(x, cov) {
	this.lm <- lm(x ~ cov);
	this.anova <- anova(this.lm);
	this.result <- matrix(NA, nrow = 1, ncol = 4);
	colnames(this.result) <- c("sum.of.squares", "mean.square", "f.value", "p.value");
	this.result[1,1] <- this.anova$`Sum Sq`[1];
	this.result[1,2] <- this.anova$`Mean Sq`[1];
	this.result[1,3] <- this.anova$`F value`[1];
	this.result[1,4] <- this.anova$`Pr(>F)`[1];
	return(this.result);
}

gse <- getGEO('GSE40005', destdir=".", GSEMatrix=TRUE);
conds <- array(c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1), dim=c(24));

anova <- apply(na.omit(exprs(gse[[1]])), 1, anova.matrix, cov = factor(conds));
anova <- data.frame(t(anova));
colnames(anova) <- c("sum.of.squares", "mean.square", "f.value", "p.value");
# Store calculations into anova data frame
anova$p.value.adj <- p.adjust(anova$p.value, method = "BH");
# Add id column, change rownames to numbers, reorder columns
anova$id <- rownames(anova);
rownames(anova) <- 1:nrow(anova);
anova <- anova[,c(6,1:5)];
sort.anova <- anova[order(anova$p.value.adj),]; # sort results by p-values

pval_row_index <- which((anova$p.value.adj < 0.05) %in% c(TRUE));
filtered_eset <- exprs(gse[[1]])[pval_row_index, ]; # get normalized data with p.value.adj < 0.05
# mean_row_index <- which((rowMeans(filtered_norm_rna) > ???) %in% c(TRUE));
# filtered_norm_rna <- filtered_norm_rna[mean_row_index, ]; # get normalized data with mean > ???

# Convert _ to whitespace
col_names <- colnames(data.matrix(filtered_eset));
col_names <- gsub("_", " ", col_names);
colnames(filtered_eset) <- col_names;
# Adjust the margins
row_names <- rownames(data.matrix(filtered_eset));
max_char_row <- max(nchar(row_names));
max_char_col <- max(nchar(col_names));
# Adjust image size
# size for plot + size for label margins + size for dendrogram in pixels + outer margins
heatmap_height <- length(row_names) * 30 + max_char_col * 8 + 50 + 50;
heatmap_width <- length(col_names) * 30 + max_char_row * 8 + 50 + 50 + 75; # add 75 for legend
png("heatmap.png", height = heatmap_height, width = heatmap_width, res = 72);
pheatmap(filtered_eset, cluster_cols = TRUE, legend = FALSE, scale = "row", border_color = NA, cellwidth = 30, cellheight = 30, treeheight_row = 75, treeheight_col = 75, annotation_legend = FALSE, fontsize = 12, fontfamily = "mono", fontface = "plain", color = colorRampPalette(c("green", "black", "red"))(100));
dev.off();