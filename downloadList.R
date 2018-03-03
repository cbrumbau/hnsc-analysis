# Load required packages
library("GEOquery"); # for getGEO

data = readLines("geolist.txt")
for (geo in data){
	if (nchar(geo) > 0) {
		getGEOfile(geo, destdir='.', AnnotGPL=TRUE, amount=c("full"))
	}
}