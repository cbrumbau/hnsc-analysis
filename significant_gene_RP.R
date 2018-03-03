library(GEOquery)
library(RankProd)

results <- list.files(path=".", pattern="*.RData");
for (RData in results) {
	load(RData);
	geo <- strsplit(RData, "\\.")[[1]][1]
	this.geo <- getGEO(strsplit(geo, "\\.")[[1]][1], destdir=".", GSEMatrix=TRUE);
	# Check if single case or more
	#if (ncol(RP.out$pfp) == 1) {
	#	# Single case
	#}
	# Get statisticially significant genes by PFP cutoff 0.05
	results <- topGene(RP.out, cutoff=0.05, method="pfp");
	# Map array to gene symbols
	probe_ids <- rownames(results$Table2);
	if ("Gene Symbol" %in% colnames(fData(this.geo[[1]]))) {
		gene_symbol <- fData(this.geo[[1]])[,c("ID", "Gene Symbol")];
		converted_results <- as.character(gene_symbol$"Gene Symbol"[match(probe_ids, gene_symbol$"ID")]);
	}
	if ("GENE_SYMBOL" %in% colnames(fData(this.geo[[1]]))) {
		gene_symbol <- fData(this.geo[[1]])[,c("ID", "GENE_SYMBOL")];
		converted_results <- as.character(gene_symbol$"GENE_SYMBOL"[match(probe_ids, gene_symbol$"ID")]);
	}
	if ("miRNA_ID" %in% colnames(fData(this.geo[[1]]))) {
		mirna_id <- fData(this.geo[[1]])[,c("ID", "miRNA_ID")];
		converted_results <- as.character(mirna_id$"miRNA_ID"[match(probe_ids, mirna_id$"ID")]);
	}
	if ("gene_assignment" %in% colnames(fData(this.geo[[1]]))) {
		gene_assignment <- fData(this.geo[[1]])[,c("ID", "gene_assignment")];
		converted_results <- as.character(gene_assignment$"gene_assignment"[match(probe_ids, gene_assignment$"ID")]);
	}
	if ("GB_RANGE" %in% colnames(fData(this.geo[[1]]))) {
		gb_range <- fData(this.geo[[1]])[,c("ID", "GB_RANGE")];
		converted_results <- as.character(gb_range$"GB_RANGE"[match(probe_ids, gb_range$"ID")]);
	}
	# Store gene symbols to files
	try(write(converted_results, file=paste(geo,".txt", sep=""), ncolumns=1));
	rm(converted_results);
}