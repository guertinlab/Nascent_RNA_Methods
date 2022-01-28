#!/usr/bin/env Rscript

# assumes the naming convention CELL_COND_TIME_ETC_rep<#>_<plus><minus>.bigWig
# assumes that you are running the script in the directory with your bigWigs 

require(bigWig)

coverage = c()
file.name = c()

for (plus.bigWig in Sys.glob(file.path("*_rep1_plus.bigWig"))) {
	file.prefix = strsplit(plus.bigWig, '_rep1')[[1]][1]
	#print(file.prefix)
    minus.bigWig = paste0(strsplit(plus.bigWig, 'plus')[[1]][1], "plus.bigWig")
    #print(minus.bigWig)
    }
for (replicate.plus.bigWig in Sys.glob(file.path("*_rep*_plus.bigWig"))) { 
   	#print(replicate.plus.bigWig)
   	replicate.prefix = strsplit(replicate.plus.bigWig, '_plus.bigWig')[[1]][1]
	replicate.minus.bigWig = paste0(strsplit(replicate.prefix, '_plus.bigWig')[[1]][1], "_minus.bigWig")
	#print(replicate.minus.bigWig)
	bigwig.plus = load.bigWig(replicate.plus.bigWig)
	bigwig.minus = load.bigWig(replicate.minus.bigWig)
	coverage.bigwig.plus = bigwig.plus$basesCovered*bigwig.plus$mean
	coverage.bigwig.minus = abs(bigwig.minus$basesCovered*bigwig.minus$mean)
	coverage = append(coverage, coverage.bigwig.plus)
	coverage = append(coverage, coverage.bigwig.minus)
	file.name = append(file.name, replicate.plus.bigWig)
	file.name = append(file.name, replicate.minus.bigWig)
	#print(coverage)
	unload.bigWig(bigwig.plus)
	unload.bigWig(bigwig.minus)
   	}	
  
celltype.prefix = strsplit(file.prefix, '_')[[1]][1]
	
write.table(data.frame(file.name, coverage), file = paste0(celltype.prefix, '_normalization.txt'), sep = '\t', row.names = FALSE, col.names=FALSE, quote =FALSE)

