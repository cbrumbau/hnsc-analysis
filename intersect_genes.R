# results <- list.files(path=".", pattern="*.txt");

# for (text_file in results) {
	# geo <- strsplit(RData, "\\.")[[1]][1];
	# all_lines <- scan(file=text, what=list(""), blank.lines.skip=TRUE);
# }

GSE23036 <- scan(file="GSE23036.txt", what=list(""), blank.lines.skip=TRUE);
GSE55546 <- scan(file="GSE55546.txt", what=list(""), blank.lines.skip=TRUE);
GSE55547 <- scan(file="GSE55547.txt", what=list(""), blank.lines.skip=TRUE);
GSE55548 <- scan(file="GSE55548.txt", what=list(""), blank.lines.skip=TRUE);
GSE55549 <- scan(file="GSE55549.txt", what=list(""), blank.lines.skip=TRUE);

gene.list <- list(GSE23036, GSE55546, GSE55547, GSE55548, GSE55549)
Reduce("intersect", gene.list)