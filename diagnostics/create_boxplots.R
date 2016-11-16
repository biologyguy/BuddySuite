
setwd("/Users/depereirabarrsn/Buddysuite/diagnostics/results")


### Read in seqbuddy avg times
file <- "sb_avg_all_times.txt"
sb <- as.matrix(read.table(file, header=TRUE, sep = "\t", row.names = 1, as.is=TRUE))

file <- "ab_avg_all_times.txt"
ab <- as.matrix(read.table(file, header=TRUE, sep = "\t", row.names = 1, as.is=TRUE))

file <- "pb_avg_all_times.txt"
pb <- as.matrix(read.table(file, header=TRUE, sep = "\t", row.names = 1, as.is=TRUE))

n <- c(" ", " 10 ", " 100 ", " 1000 ", " 10000", " ")
azul <- rgb(114,147,203,255, maxColorValue=255)
verde <- rgb(132,186,91,255, maxColorValue=255)
vermelho <- rgb(211,94,96,255, maxColorValue=255)
black <- rgb(255,255,255,255, maxColorValue=255)

ats <- c(0, 1, 10, 100, 1000)

p1 <- pb[,1]
a1 <- ab[,1]
s1 <- sb[,1]
p2 <- pb[,2]
a2 <- ab[,2]
s2 <- sb[,2]
p3 <- pb[,3]
a3 <- ab[,3]
s3 <- sb[,3]
p4 <- pb[,4]
a4 <- ab[,4]
s4 <- sb[,4]

z <- c(NA, NA) ## for spacing the columns

length(z) = length(p1) = length(a1) = length(s1)
length(p2) = length(a2) = length(s2)
length(p3) = length(a3) = length(s3)
length(p4) = length(a4) = length(s4)

cl <- c(azul, black, black, black, verde, black, black, black, vermelho, black, black, black, black, black, black, azul, black, black, black, verde, black, black, black, vermelho, black, black, black, black, black, black, azul, black, black, black, verde, black, black, black, vermelho, black, black, black, black, black, black, azul, black, black, black, verde, black, black, black, vermelho)

tempos <- data.frame(s1, z, z, z, a1, z, z, z, p1, z, z, z, z, z, z, s2, z, z, z, a2, z, z, z, p2, z, z, z, z, z, z, s3, z, z, z, a3, z, z, z, p3, z, z, z, z, z, z, s4, z, z, z, a4, z, z, z, p4, z)

svg("/Users/depereirabarrsn/BuddySuite/diagnostics/results/figure.svg")
## bottom, left, top, right
par(mar=c(5,6,2,3)+0.1)
#boxplot(sb[,1], ab[,1], pb[,1], sb[,2], ab[,2], pb[,2], sb[,3], ab[,3], pb[,3], sb[,4], ab[,4], pb[,4], at = c(1,2,3,6,7,8,11,12,13,16,17,18), col= cl, axes = FALSE, log = "y", border = cl)
stripchart(tempos, method="jitter", jitter=1.2,  vertical=TRUE, col= cl, axes = FALSE, log = "y", pch = 16)
axis(1, at=c(0,5,20,35,50,56), labels=n)
axis(2, at=ats, labels=paste(ats, "", sep=""), las =2)
mtext(side = 2, "Average time in seconds", line = 3.5)
mtext(side = 1, "Number of sequences", line = 3)
legend(1, 3000, legend = c("SeqBuddy", "AlignBuddy", "PhyloBuddy"), fill = c(azul, verde, vermelho), border = c(azul, verde, vermelho), bty = "o")
dev.off()





