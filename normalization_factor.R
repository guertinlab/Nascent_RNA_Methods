#!/usr/bin/env Rscript

# assumes the naming convention CELL_COND_TIME_ETC_rep<#>_<plus><minus>.bigWig for unpaired reads or
# CELL_COND_TIME_ETC_rep<#>_<plus><minus>_PE<1><2>.bigWig for paired reads
# assumes that you are running the script in the directory with your bigWigs 



require(bigWig)

coverage = c()
file.name = c()

# Establishes cell type for naming of output (see end of script)
for (plus.bigWig in Sys.glob(file.path("*_rep1_plus*bigWig"))) {
	file.prefix = strsplit(plus.bigWig, '_rep1')[[1]][1]
	#print(file.prefix)
    #minus.bigWig = paste0(strsplit(plus.bigWig, 'plus')[[1]][1], "plus.bigWig") 
    #print(minus.bigWig)
    }

# Creates one coverage dataframe for all samples in working directory
for (replicate.plus.bigWig in Sys.glob(file.path("*_rep*_plus*bigWig"))) { 
   	#print(paste("Current replicate:",replicate.plus.bigWig))
   	
   	#get suffix as variable 
   	replicate.suffix = strsplit(replicate.plus.bigWig, "_plus")[[1]][2]
   	replicate.pair = strsplit(replicate.suffix, ".bigWig")[[1]][1]
   	minus.suffix = paste0("_minus", replicate.pair, ".bigWig")
   	# Get name for minus files
   	replicate.prefix = strsplit(replicate.plus.bigWig, '_plus')[[1]][1]
	replicate.minus.bigWig = paste0(replicate.prefix, minus.suffix)
	#print(replicate.minus.bigWig)
	
	# Load bigWigs
	bigwig.plus = load.bigWig(replicate.plus.bigWig)
	bigwig.minus = load.bigWig(replicate.minus.bigWig)
	
	# Calculate coverage
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


# Name according to cell type and export as txt file
celltype.prefix = strsplit(file.prefix, '_')[[1]][1]	
write.table(data.frame(file.name, coverage), file = paste0(celltype.prefix, '_normalization.txt'), sep = '\t', row.names = FALSE, col.names=FALSE, quote =FALSE)
