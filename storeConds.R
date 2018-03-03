all.names <- scan("geolist.txt", what="", sep="\n")
all.conds <- list(
array(c(rep(0,73),rep(1,23)), dim=c(96)),
array(c(rep(c(0,1),22)), dim=c(44)),
array(c(rep(0,28),rep(1,8)), dim=c(36)),
array(c(rep(1,270)), dim=c(270)),
array(c(rep(1,134)), dim=c(134)),
array(c(rep(1,8)), dim=c(8)),
array(c(rep(c(1,0),5)), dim=c(10)),
array(c(rep(c(0,1),15)), dim=c(30)),
array(c(rep(0:1,each=3)), dim=c(6)),
array(c(rep(1,13),rep(0,5)), dim=c(18)),
array(c(rep(0,4),rep(1,8)), dim=c(12)),
array(c(rep(1,4),rep(0,4)), dim=c(8)),
array(c(rep(1,16),rep(0,4)), dim=c(20)),
array(c(rep(1,12),rep(0,4)), dim=c(16)),
array(c(rep(1,24)), dim=c(24)),
array(c(rep(1,16)), dim=c(16)),
array(c(rep(1,225)), dim=c(225)),
array(c(rep(1,138)), dim=c(138)),
array(c(rep(1,44),rep(0,25)), dim=c(69)),
array(c(rep(1,44),rep(0,25)), dim=c(69)),
array(c(rep(1,44),rep(0,25)), dim=c(69)),
array(c(rep(1,6)), dim=c(6)),
array(c(rep(1,42)), dim=c(6)),
array(c(rep(1,6)), dim=c(6)),
array(c(rep(1,316)), dim=c(316)),
array(c(1,1,0,0), dim=c(4)),
array(c(rep(1,63),rep(0,5)), dim=c(68)),
array(c(rep(1,94)), dim=c(94)),
array(c(rep(1,10)), dim=c(10)),
array(c(rep(1,24)), dim=c(24)),
array(c(rep(1,44)), dim=c(44)),
array(c(rep(1,39)), dim=c(39)),
array(c(rep(1:0,each=16)), dim=c(32))
)
conds.df <- data.frame(all.names, row.names=all.names)
conds.df$conds <- all.conds
conds.df <- subset(conds.df, select=c('conds'))