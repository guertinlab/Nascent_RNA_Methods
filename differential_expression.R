#!/usr/bin/env Rscript

Args=commandArgs(TRUE)
require(lattice)
require(DESeq2)
n.untreated = 2
n.treated =2


plot.ma.lattice <- function(ma.df, filename = 'file.name', p = 0.01, title.main = "Differential PRO Expression",log2fold = 0.5)
  {
  pdf(paste("MA_plot_", filename, ".pdf", sep=''), width=3.83, height=3.83)
  print(xyplot(ma.df$log2FoldChange ~ log(ma.df$baseMean, base=10),
               groups=(ma.df$padj < p & abs(ma.df$log2FoldChange) > log2fold & !is.na(ma.df$padj)),
               col=c("grey40","red"), main=title.main, scales="free", aspect=1, pch=20, cex=0.5,
               ylab=expression("log"[2]~"PRO fold change"), xlab=expression("log"[10]~"Mean of Normalized Counts"),
               par.settings=list(par.xlab.text=list(cex=1.1,font=2), par.ylab.text=list(cex=1.1,font=2))))
  dev.off()
}




x = read.table(Args[1], sep = '\t', header = TRUE)
rownames(x) = x[,1]

x = x[,seq(2,to=ncol(x),by=2)]

x = x[,c(grep(Args[2], colnames(x)), c(1:length(colnames(x)))[!(c(1:length(colnames(x))) %in% grep(Args[2], colnames(x)))])]

sample.conditions = factor(c(rep("untreated",n.untreated), rep("treated",n.treated)),
                           levels=c("untreated","treated"))

deseq.df = DESeqDataSetFromMatrix(x, DataFrame(sample.conditions), ~ sample.conditions)

deseq.df =estimateSizeFactors(deseq.df)
deseq.df = estimateDispersions(deseq.df)
deseq.df = nbinomWaldTest(deseq.df)
res.deseq = results(deseq.df)
plot.ma.lattice(res.deseq, filename = Args[3], p = 0.05,log2fold = 0.0)
