#!/usr/bin/env Rscript

Args=commandArgs(TRUE)

exon = read.table(Args[1], header =FALSE)
intron = read.table(Args[2], header =FALSE)

prefix.name = strsplit(Args[1], "_exon_counts.bed")[[1]][1]


exon.data = data.frame(aggregate(V8 ~ V4, data = exon, sum), aggregate(V7 ~ V4, data = exon, sum))[,c(1,2,4)]
colnames(exon.data) = c('gene', 'exon.length', 'exon.signal')

intron.data = data.frame(aggregate(V8 ~ V4, data = intron, sum), aggregate(V7 ~ V4, data = intron, sum))[,c(1,2,4)]
colnames(intron.data) = c('gene', 'intron.length', 'intron.signal')

all.data = merge(intron.data, exon.data, by = 'gene')

all.data$exon.density = all.data$exon.signal/all.data$exon.length
all.data$intron.density = all.data$intron.signal/all.data$intron.length

all.data$ratio = all.data$exon.density / all.data$intron.density

pdf(paste(prefix.name,'_exon_intron_ratio.pdf', sep=''), width=2.5, height=4, useDingbats=FALSE)
par(pty="s")
boxplot(log(all.data$ratio, base = 10), main = paste('median = ', round(median(all.data$ratio, na.rm = TRUE), 2)),
        pch = 16, outline=FALSE, ylab=expression("log"[10]~"(Exon Density / Intron Density)"))
abline(h = log(1, base = 10), col='blue', lty = 2)
dev.off()

