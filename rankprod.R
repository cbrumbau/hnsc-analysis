# Load required packages
library("GEOquery"); # for getGEO
library("RankProd"); # for RPadvance

gse <- getGEO('GSE40005', destdir=".", GSEMatrix=TRUE);
conds <- array(c(0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1), dim=c(24));
origin <- rep(1,24);

RP.out <- RPadvance(exprs(gse[[1]]), conds, origin, num.perm=100, na.rm=FALSE, plot=FALSE, rand=123, huge=TRUE);