# Load required packages
library("GEOquery"); # for getGEO
library("RankProd"); # for RPadvance

source("storeConds.R")
data = readLines("geolist.txt")
for (geo in data) {
	# Read in data
	if (nchar(geo) > 0) {
		if (grepl("GSE", geo)) {
			this.geo <- getGEO(geo, destdir=".", GSEMatrix=TRUE);
		} else {
			this.geo <- getGEO(geo, destdir=".");
		}
		# Store test conditions
		this.conds <- conds.df[geo,][[1]]
		this.origin <- rep(1, length(this.conds))
		# Run, store data and reset
		if (grepl("GSE", geo)) {
			try(RP.out <- RPadvance(exprs(this.geo[[1]]), this.conds, this.origin, num.perm=100, na.rm=FALSE, plot=FALSE, rand=123, huge=TRUE));
		} else {
			try(RP.out <- RPadvance(Table(this.geo), this.conds, this.origin, num.perm=100, na.rm=FALSE, plot=FALSE, rand=123, huge=TRUE));
		}
		file.name <- paste(geo, ".RData", sep='');
		try(save(RP.out, file=file.name));
		try(rm(RP.out));
		rm(this.geo);
	}
}