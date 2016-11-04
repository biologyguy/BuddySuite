#!/usr/bin/Rscript

setwd("/Users/depereirabarrsn/BuddySuite/diagnostics")

### Read in seqbuddy avg times
file <- "sb_avg_times.txt"
sb <- as.matrix(read.table(file, header=TRUE, sep = "\t", row.names = 1, as.is=TRUE))

file <- "ab_avg_times.txt"
ab <- as.matrix(read.table(file, header=TRUE, sep = "\t", row.names = 1, as.is=TRUE))

file <- "pb_avg_times.txt"
pb <- as.matrix(read.table(file, header=TRUE, sep = "\t", row.names = 1, as.is=TRUE))

n <- c("", " 10 ", " 100 ", " 1000 ", " ")
azul <- rgb(114,147,203,255, maxColorValue=255)
verde <- rgb(132,186,91,255, maxColorValue=255)
vermelho <- rgb(211,94,96,255, maxColorValue=255)
cl <- c(azul, verde, vermelho, azul, verde, vermelho, azul, verde, vermelho)
# cl <- c("deepskyblue3","palegreen3", "deepskyblue3","palegreen3", "deepskyblue3","palegreen3")

ats <- c(0, 1, 10, 100, 1000)

#svg("temp.svg")
## bottom, left, top, right
par(mar=c(5,6,2,3)+0.1)
# boxplot(log10(sb[,1]), log10(ab[,1]), log10(pb[,1]), log10(sb[,2]), log10(ab[,2]), log10(pb[,2]), log10(sb[,3]), log10(ab[,3]), log10(pb[,3]), at = c(1,2,3,5,6,7,9,10,11), col= cl, axes = FALSE, log = "y")
boxplot(sb[,1], ab[,1], pb[,1], sb[,2], ab[,2], pb[,2], sb[,3], ab[,3], pb[,3], at = c(1,2,3,7,8,9,13,14,15), col= cl, axes = FALSE, log = "y", border = cl)
axis(1, at=c(0,2,8,14,16), labels=n)
axis(2, at=ats, labels=paste(ats, "", sep=""), las =2)
mtext(side = 2, "Average time in seconds log1010 scale", line = 3.5)
mtext(side = 1, "Number of sequences", line = 3)
legend(1, 2000, legend = c("SeqBuddy", "AlignBuddy", "PhyloBuddy"), fill = c(azul, verde, vermelho), border = c(azul, verde, vermelho), bty = "o")
#dev.off()





