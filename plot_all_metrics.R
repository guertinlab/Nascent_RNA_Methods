#!/usr/bin/env Rscript
Args=commandArgs(TRUE)

require(lattice)
require(grid)

x = read.table(Args[1], sep = '\t', header = TRUE)

if (length(unique(x$experiment)) == 1) {
    width.pdf = 2.3
} else {
    width.pdf = 0.8219444+(0.7945139*length(unique(x$experiment)))
    }

pdf(paste0(Args[2], '.pdf'), width = width.pdf, height=10)
  print(
      barchart
      (value ~ experiment | metric, data = x,
          par.settings=list(par.xlab.text=list(cex=0.9,font=1),
                            par.ylab.text=list(cex=0.7,font=1),
               strip.background=list(col="grey80", cex = 0.6),
               par.main.text=list(cex=0.85, font=1)),
          layout=c(1,7),
          cex.axis=0.1,
               par.strip.text=list(cex=0.8, font=1, col='black'),
         ylim = list(
             c(0,1),
             c(0,1),
             c(0,max(c(x$value[x$metric == 'Degradation Ratio'],x$threshold[x$metric == 'Degradation Ratio'])) + max(c(x$value[x$metric == 'Degradation Ratio'],x$threshold[x$metric == 'Degradation Ratio']))/10),
             c(0,max(c(x$value[x$metric == 'Exon/Intron Ratio'],x$threshold[x$metric == 'Exon/Intron Ratio'])) + max(c(x$value[x$metric == 'Exon/Intron Ratio'],x$threshold[x$metric == 'Exon/Intron Ratio']))/10),
             c(0,max(c(x$value[x$metric == 'Pause Index'],x$threshold[x$metric == 'Pause Index'])) + max(c(x$value[x$metric == 'Pause Index'],x$threshold[x$metric == 'Pause Index']))/10),
             c(0,max(c(x$value[x$metric == 'rDNA Alignment Rate'],x$threshold[x$metric == 'rDNA Alignment Rate'])) + max(c(x$value[x$metric == 'rDNA Alignment Rate'],x$threshold[x$metric == 'rDNA Alignment Rate']))/10),
             c(0,1)),
                   col='grey90',between=list(y=0.5, x=0.5),
         ylab="    ",
         xlab = 'Experiment',
         scales=list(x=list(cex=0.8, rot = 30), y =list(cex=0.8, relation="free", axs = "i")),
         threshold = x$threshold,
                   panel=function(..., threshold, subscripts, fill.color, fill.color.2) {
                       threshold = unique(threshold[subscripts])
                       fill.color = c('#D2F9D2','#FF5076', '#D2F9D2', '#D2F9D2', '#FF5076', '#D2F9D2', '#FF5076')
                       fill.color.2 = c('#FF5076','#D2F9D2', '#FF5076', '#FF5076', '#D2F9D2', '#FF5076', '#D2F9D2')
                       grid.rect(x=unit(0.5, "npc"), y=unit(threshold, "native"),
                                 width=unit(1, "npc"), height=unit(1, "npc"),
                                 just = c("center", "top"),
                                 gp = gpar(col = 'transparent', lty=2, fill = fill.color[packet.number()]))
                       grid.rect(x=unit(0.5, "npc"), y=unit(threshold, "native"),
                                 width=unit(1, "npc"), height=unit(1, "npc"),
                                 just = c("center", "bottom"),
                                 gp = gpar(col = 'transparent', lty=2, fill = fill.color.2[packet.number()]))
                       panel.barchart(...)
                       panel.abline(h=threshold, lty= 2, col = 'grey20')
                   }
         )
  )
dev.off()



