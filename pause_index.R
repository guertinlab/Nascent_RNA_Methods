#!/usr/bin/env Rscript

Args=commandArgs(TRUE)
argsLen <- length(Args)

densities = read.table(Args[1], header =FALSE)
prefix.name = strsplit(Args[1], "_pause_body.bed")[[1]][1]
median(densities[,10] )

pdf(paste(prefix.name,'_pause_index.pdf', sep=''), width=2.5, height=4, useDingbats=FALSE)
par(pty="s")
boxplot(log(densities[,10], base = 10), main = paste('median = ', round(median(densities[,10] ), 2)),
        pch = 16, outline=FALSE, ylab=expression("log"[10]~"(Pause Index)"))
abline(h = log(10, base = 10), col='blue', lty = 2)
dev.off()

