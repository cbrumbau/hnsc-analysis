FC.cut <- 4;
FDR.cut <- 10^-25;
up.gene <- c();
down.gene <- c();
for (idx in 1:nrow(data_diffexp.tab)) {
	if (data_diffexp.tab$FDR[idx] < FDR.cut) {
		if (data_diffexp.tab$logFC[idx] > FC.cut) {
			up.gene <- c(up.gene, rownames(data_diffexp.tab)[idx]);
		}
		if (data_diffexp.tab$logFC[idx] < FC.cut) {
			down.gene <- c(down.gene, rownames(data_diffexp.tab)[idx]);
		}
	}
}
write(up.gene, file = "upregulated_genes.txt", append = "FALSE");
write(down.gene, file = "downregulated_genes.txt", append = "FALSE");