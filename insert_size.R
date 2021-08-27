#!/usr/bin/env Rscript

Args=commandArgs(TRUE)
insert.histogram = Args[1]
#ratio.degraded.intact = Args[2]
umi.length = as.numeric(Args[2])


prefix.name = strsplit(insert.histogram, ".hist")[[1]][1]
hist = read.table(insert.histogram, header=FALSE)
colnames(hist) = c('size', 'frequency')
hist$size = hist$size - umi.length

degraded = sum(hist$frequency[(10 <= hist$size) & (hist$size <= 20)])
intact = sum(hist$frequency[(30 <= hist$size) & (hist$size <= 40)])

ratio.degraded.intact = degraded/intact
  

ratio.degraded.intact.2 = round(as.numeric(ratio.degraded.intact), 2)

system(paate0('echo -e "',ratio.degraded.intact.2,'\t',prefix.name,'\t1.0\tDegradation Ratio" >> prefix.name_QC_metrics.txt'))

pdf(paste(prefix.name, '_sizes', '.pdf', sep=''), width=5, height=5, useDingbats=FALSE)
par(pty="s")
plot(frequency~size, hist, 
     pch = 16, cex = 0.5, ylab = 'Frequency', xlim = c(0, max(hist$size)), xlab = 'Insert Size', xaxs="i",
     main = bquote('Ratio of degraded to intact RNA ='~.(ratio.degraded.intact.2)))
abline(v = 10, col='blue', lty = 2)
abline(v = 20, col='blue', lty = 2)
abline(v = 30, col='red', lty = 2)
abline(v = 40, col='red', lty = 2)
#text(max(hist$size), max(hist$frequency)*0.95, bquote('Ratio of degraded to intact RNA = '~.(ratio.degraded.intact.2)),  pos = 2)

dev.off()

