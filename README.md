# Nascent_RNA_Methods

Precision genomic run-on assays (PRO-seq) quantify nascent RNA at single nucleotide resolution with strand specificity. Here, we deconstruct a recently published genomic nascent RNA processing pipeline (PEPPRO) into its components. The analyses are presented as individual code chunks with comprehensive details. PRO-seq experiments are evolving and variations can be found throughout the literature. We describe analysis of paired end libraries that contain a unique molecular identifier on the the 5´ adapter adjacent to the RNA insert. These two features allow for calculations of RNA fragment length and library complexity. Our hope is that the end users will adopt the analysis framework and modify the workflow as needed. We describe and provide the code to calculate the following quality control metrics: library complexity, nascent RNA purity, nuclear run on efficiency, alignment rate, sequencing depth, and RNA integrity.        

# One time [better word for one-time] software installations, downloads, and processing steps

## Dependencies, software, and scripts

Although we present novel quality control metrics and specialized software herein, most of the workflow depends upon more general software. Fortunately, this software is well-maintained and documented, so we only provide a short description and the links below. 

[Cite all the software in line]

`bedtools` :  a comprehensive suite of tools that efficiently perform a wide range of operations on genomic intervals.  https://bedtools.readthedocs.io/en/latest/

`bowtie2` : aligns sequence reads to reference sequences. http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

`cutadapt` : removes a defined sequence, such as adapter sequence, from sequence reads.  https://cutadapt.readthedocs.io/en/stable/

`fastq_pair` : outputs only sequence reads that have a matched paired end read and a separates unpaired reads. https://github.com/linsalrob/fastq-pair 

`FLASH` : merges paired end reads by detectign overlap. https://ccb.jhu.edu/software/FLASH/ 

`fqdedup` : remove duplicated sequences from FASTQ files. https://github.com/guertinlab/fqdedup

`samtools` : a suite of tools for parsing and interfacing with high throughput sequencing data files. http://www.htslib.org

`seqOutBias` : software designed to correct enzymatic sequence biases that also has options to output BED and bigWig files with desired features. https://github.com/guertinlab/seqOutBias

`seqtk` : a multifunctional toolkit for processing sequence files. We only use two functionalities: trimming a defined number of bases from the ends of reads and reverse complementing sequence reads. https://github.com/lh3/seqtk

`wget` : retrieves files from a wide range of internet protocols. https://www.gnu.org/software/wget/

In addition, to facilitate data analysis and graphical output, we developed the software and R scripts. The following code uses `wget` to retrieve the software and scripts. The command `chmod +x` changes the permissions of the files to executable. 


```
wget https://raw.githubusercontent.com/guertinlab/fqComplexity/main/fqComplexity
wget https://raw.githubusercontent.com/guertinlab/fqComplexity/main/complexity_pro.R
wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/insert_size.R
wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/pause_index.R
wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/exon_intron_ratio.R

chmod +x insert_size.R
chmod +x fqComplexity
chmod +x complexity_pro.R
chmod +x pause_index.R
chmod +x exon_intron_ratio.R
```
The software dependencies and R scripts should be moved to a directory within the `$PATH` variable.      

## Reference genomes

PRO-seq experiments have been performed in a variety of organisms including yeast (cite Booth), Drosophila (cite Kwak and cite Duarte), and humans (cite Sathyan). Analysis of the data requires alignment to a reference genome annotation. The first step is to use `wget` to retrieve the zipped reference genome. Many websites host the assembly data and one can simply use Google the organism name combined with the term "reference genome" to find the appropriate FASTA file. The `gunzip` command unzips the refence genome file and `bowtie2-build` indexes the file to allow for efficient alignment. The code also retrieves, unzips, and builds the human rDNA reference genome (cite refgenie). 

```
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
bowtie2-build hg38.fa hg38

wget https://github.com/databio/ref_decoy/raw/master/human_rDNA.fa.gz
gunzip human_rDNA.fa.gz
bowtie2-build human_rDNA.fa human_rDNA
```

## Reference gene annotation

The quality control metrics outlined herein require the counting of sequence reads that align to three genomic features: exons, intron, and promoter-proximal pause regions. Gene annotations are available from many sources and we outline retrieval and parsing of GTF files from Ensembl (cite). The Ensembl website (http://www.ensembl.org/index.html) contains the information for the latest release, which at the time of writing this manuscript is release 104 for hg38. After retrieving and unzipping the file we parse out all exon 1 annotations--note that a single gene can have multiple first exons due to the presence of different gene isoforms. Ensembl chromosome numbers do not include the preceding "chr", so the first `sed` command appends "chr" to the chromosome name. The output of this is piped to `awk`, which prints the indicated fields. Subsequent `sed` commands drop the semicolon and quote characters from the gene and Ensembl IDs while editing the mitochondrial chromosome to match the reference genome, "chrM" as opposed to "chrMT". Finally the exon output is sorted by the first, then second column, in ascending order. The gene annotations are sorted by gene name. The resultant BED6 files for the exons include the chromosome coordinates in columns 1-3, Ensembl transcript or gene ID (ENST/ENSG), gene name, and strand information.       

NOTE: the Ensembl transcript ID ENST and Ensembl Gene ID ENSG IDs were in different columns. $14 and $10 were switched in the gene annotations, so subsequent code may break. I changed the join command to specify column 5, but I have not tested it yet.

```

release=104

wget http://ftp.ensembl.org/pub/release-${release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${release}.chr.gtf.gz
gunzip Homo_sapiens.GRCh38.${release}.chr.gtf.gz

# extract all exon 1 annotations
grep 'exon_number "1"' Homo_sapiens.GRCh38.${release}.chr.gtf | \
    sed 's/^/chr/' | \
    awk '{OFS="\t";} {print $1,$4,$5,$14,$20,$7}' | \
    sed 's/";//g' | \
    sed 's/"//g' | sed 's/chrMT/chrM/g' | sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.${release}.tss.bed

#extract all exons
grep 'exon_number' Homo_sapiens.GRCh38.${release}.chr.gtf | \
    sed 's/^/chr/' | \
    awk '{OFS="\t";} {print $1,$4,$5,$14,$20,$7}' | \
    sed 's/";//g' | \
    sed 's/"//g' | sed 's/chrMT/chrM/g' | sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.${release}.all.exons.bed

#extract all complete gene annotations
awk '$3 == "gene"' Homo_sapiens.GRCh38.${release}.chr.gtf | \
    sed 's/^/chr/' | \
    awk '{OFS="\t";} {print $1,$4,$5,$10,$14,$7}' | \
    sed 's/";//g' | \
    sed 's/"//g' | sed 's/chrMT/chrM/g' | sort -k5,5 > Homo_sapiens.GRCh38.${release}.bed
```

The goal of the following operations is to define a set of exons that excludes all instances of first exons, define all introns, and define a set of all potential pause regions for a gene by taking the region for 20 - 120 downstream of all exon 1 annotations. The `mergeBed` command collapses all overlapping intervals and the gene name information is lost. We exclude all first exon coordinates with `subtractBed`, then `intersectBed` reassigns gene names to all remaining exons. The final `awk` command defines a 100 base pause region window downstream of all transcription start sites based on the gene strand.    

```
#merge exon intervals that overlap each other
mergeBed -s -c 6 -o distinct -i Homo_sapiens.GRCh38.${release}.all.exons.bed | \
    awk '{OFS="\t";} {print $1,$2,$3,$4,$2,$4}' | 
    sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.${release}.all.exons.merged.bed

#remove all first exons (so pause region is excluded from exon / intron density ratio)
subtractBed -s -a Homo_sapiens.GRCh38.${release}.all.exons.merged.bed -b Homo_sapiens.GRCh38.${release}.tss.bed | \
    sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.${release}.no.first.exons.bed

#define and name all introns 
subtractBed -s -a Homo_sapiens.GRCh38.${release}.bed -b Homo_sapiens.GRCh38.${release}.all.exons.merged.bed | \
    sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.${release}.introns.bed 

#extract gene names of exons
intersectBed -s -wb -a Homo_sapiens.GRCh38.${release}.no.first.exons.bed -b Homo_sapiens.GRCh38.${release}.bed | \
    awk '{OFS="\t";} {print $1,$2,$3,$11,$4,$4}' | \
    sort -k1,1 -k2,2n >  Homo_sapiens.GRCh38.${release}.no.first.exons.named.bed

#extract the pause region from the first exons, position 20 - 120 downstream of the TSS
awk  '{OFS="\t";} $6 == "+" {print $1,$2+20,$2 + 120,$4,$5,$6} \
    $6 == "-" {print $1,$3 - 120,$3 - 20,$4,$5,$6}' Homo_sapiens.GRCh38.${release}.tss.bed  | \
    sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.${release}.pause.bed
    
```

Lastly, we retrieve, sort, and harmonize chromosome size information.


```
wget https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes

#sort chromosome sizes and harmonize chrM name
sort -k1,1 -k2,2n hg38.chrom.sizes | \
    sed 's/chrMT/chrM/g' > hg38.chrom.order.txt
    
```

# Processing PRO-seq data

## Initialize variables

```
#Ensembl release used above to parse gene annotations
release=104
#Length of UMI on 5' adapter
UMI_length=8
#PE1 insert size (read length minus UMI length)
read_size=30
#Number of cores to be used in bowtie2 for parallel processing
cores=6
#Directory where the files are
directory=/Users/guertinlab/Downloads/Batch1 

#File prefix if you run this in a loop
name=$(echo $1 | awk -F"_PE1.fastq.gz" '{print $1}')

#File prefix for an individual sample (remove ".fastq.gz" and PE1 designation from the file name)
name=LNCaP_10uMEnza_rep3_batch2
```

## Begin
```
cd $directory 
echo $name
gunzip ${name}_PE*.fastq.gz
```

Remove adapter sequences and inserts less than 10 bases

(We can also do --cores=$cores or just not do it in parallel if it's going to cause issues)
```
cutadapt --cores=$cores -m $((UMI_length+10)) -O 1 -a TGGAATTCTCGGGTGCCAAGG ${name}_PE1.fastq -o ${name}_PE1_noadap.fastq --too-short-output ${name}_PE1_short.fastq > ${name}_PE1_cutadapt.txt
cutadapt --cores=$cores -m $((UMI_length+10)) -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC ${name}_PE2.fastq -o ${name}_PE2_noadap.fastq --too-short-output ${name}_PE2_short.fastq > ${name}_PE2_cutadapt.txt
```

Remove PCR duplicates from PE1
```
fqdedup -i ${name}_PE1_noadap.fastq -o ${name}_PE1_dedup.fastq
```

## DEGRADATION RNA INTEGRITY
```
reads=$(wc -l ${name}_PE1_noadap.fastq | awk '{print $1/4}')
fastq_pair -t $reads ${name}_PE1_noadap.fastq ${name}_PE2_noadap.fastq

flash -q --compress-prog=gzip --suffix=gz ${name}_PE1_noadap.fastq.paired.fq ${name}_PE2_noadap.fastq.paired.fq -o ${name}
insert_size.R ${name}.hist ${UMI_length}

rm ${name}_PE*_noadap.fastq.paired.fq
```

## PROCESSING FOR ALIGNMENT

Trim the UMI and reverse complement
```
seqtk trimfq -b ${UMI_length} ${name}_PE1_noadap.fastq | seqtk seq -r - > ${name}_PE1_processed.fastq
seqtk trimfq -e ${UMI_length} ${name}_PE2_noadap.fastq | seqtk seq -r - > ${name}_PE2_processed.fastq
```

Remove reads aligning to rDNA
(I'm switching all flags to hex for consistency)
```
bowtie2 -p $cores -x human_rDNA -U ${name}_PE1_processed.fastq 2>${name}_bowtie2_rDNA.log | samtools sort -n - | samtools fastq -f 0x4 - > ${name}_PE1.rDNA.fastq

reads=$(wc -l ${name}_PE1.rDNA.fastq | awk '{print $1/4}')
fastq_pair -t $reads ${name}_PE1.rDNA.fastq ${name}_PE2_processed.fastq
```

Align to hg38

Should we change --maxins to something like 1000? Default is 500. I guess this doesn't matter if we're going forward with all mapped PE1 reads, but I still prefer concordantly aligned reads.
```
bowtie2 -p cores -x hg38 --rf -1 ${name}_PE1.rDNA.fastq.paired.fq -2 ${name}_PE2_processed.fastq.paired.fq 2>${name}_bowtie2_hg38.log | samtools view -b - | samtools sort - -o ${name}.bam
```

## rDNA ALIGNMENT RATE
```
PE1_prior_rDNA=$(wc -l ${name}_PE1_processed.fastq | awk '{print $1/4}')
PE1_post_rDNA=$(wc -l ${name}_PE1.rDNA.fastq | awk '{print $1/4}')
```

Should we echo these QC metrics and/or save them somewhere?

rRNA alignment rate (does not account for low genome alignment rates, so this can be artifically low if the concordant alignment rates are low)
If concodarnt alignmnet rate are low this supercedes rDNA alignment rate and multiplying by the inverse of the concordant alignment rate gives a better approximation of rDNA alignment rate
```
not_considering_overall_alignment_rate=$(echo "$(($PE1_prior_rDNA-$PE1_post_rDNA))" | awk -v myvar=$PE1_prior_rDNA '{print $1/myvar}')
```

alternatively, of the aligned reads, what fraction is rDNA:
this is what PEPPRO should do

extract concordant aligned reads from BAM
this is useful as a metric to know whether you want to sequence more, usually over 10 milion reads is good if you have 3+ replicates. 
```
concordant_pe1=$(samtools view -c -f 0x42 ${name}.bam)

overall_alignment_considered=$(echo "$(($PE1_prior_rDNA-$PE1_post_rDNA))" | awk -v myvar=$concordant_pe1 '{print $1/myvar}')
```

## FRACTION OF FILTERED READS THAT ARE MAPPABLE
this is a QC metric in itself. Concordant alignment rate is typically above 90% for good libraries
map is less stringently than concordantly mapped and most should map
```
map_pe1=$(samtools view -c -f 0x40 -F 0x4 ${name}.bam)
pre_alignment=$(wc -l ${name}_PE1.rDNA.fastq.paired.fq | awk '{print $1/4}')
alignment_rate=$(echo "scale=2 ; $map_pe1 / $pre_alignment" | bc)
```

## COMPLEXITY AND THEORETICAL READ DEPTH

calculate PE1 total raw reads
```
PE1_total=$(wc -l ${name}_PE1.fastq | awk '{print $1/4}')
```

calculate PE1 reads without adapters 
```
PE1_noadap=$(wc -l ${name}_PE1_noadap.fastq | awk '{print $1/4}')
```

This inverse of this factor is a QC metric for percent adapter/adapter ligation products (including 1 base inserts)

```
factorX=$(echo "scale=2 ; $PE1_total / $PE1_noadap" | bc)

echo fraction of reads that are not adapter/adapter ligation products
echo $factorX | awk '{print 1/$1}'
```

calculate PE1 deduplicated reads
```
PE1_dedup=$(wc -l ${name}_PE1_dedup.fastq | awk '{print $1/4}')
```

divide
```
factorY=$(echo "scale=2 ; $concordant_pe1 / $PE1_dedup" | bc)
```

two different curves:
this curve lets you know the quality of the library in terms of it's complexity.
if at 10 milion reads depth, 75% of the reads are unique, then it passes. The higher the better 

```
fqComplexity -i ${name}_PE1_noadap.fastq 
```

This curve is similar, but the goal is to estimate the raw read depth needed to acieve a target concordant aligned read count
the two factors needed are the fraction of the total reads that are adapter/adapter ligation products (this is the only step that preceeds deduplication)
 and the fraction of deduplicated reads that align concordantly after filtering rDNA-aligned reads, short reads, and unaligned reads
 simply solve for read_depth after providing a desired Concordant Aligned
rearrange the equation for them?
empirically noticed that copying and pasting the equation into R interprets the minus signs as hyphens

Change to not have it re-do the subsampling and deduplicating
```
fqComplexity -i ${name}_PE1_noadap.fastq -x factorX -y $factorY
```

counterintuitively, you can have a high quality and complex library that is not practical to sequence to further depth because the number of adapter/adapter
reads is just too high
the post-depulication factors are all individual QC metrics, so these should be considered individually to determine whether the library is high quality.


## Get the reads in a BED

I changed the flags to hex
```
samtools view -b -f 0x40 ${name}.bam > ${name}_PE1.bam
samtools view -bh -F 0x14 ${name}_PE1.bam > ${name}_PE1_plus.bam
samtools view -bh -f 0x10 ${name}_PE1.bam > ${name}_PE1_minus.bam
    
seqOutBias hg38.fa ${name}_PE1_plus.bam --no-scale --bed ${name}_PE1_plus.bed --bw=${name}_PE1_plus.bigWig --tail-edge --read-size=$read_size
seqOutBias hg38.fa ${name}_PE1_minus.bam --no-scale --bed ${name}_PE1_minus.bed --bw=${name}_PE1_minus.bigWig --tail-edge --read-size=$read_size

awk '{OFS="\t";} {print $1,$2,$3,$4,$5,"+"}' ${name}_PE1_plus.bed > ${name}_PE1_plus_strand.bed
awk '{OFS="\t";} {print $1,$2,$3,$4,$5,"-"}' ${name}_PE1_minus.bed > ${name}_PE1_minus_strand.bed

cat ${name}_PE1_plus_strand.bed ${name}_PE1_minus_strand.bed | sort -k1,1 -k2,2n > ${name}_PE1_signal.bed
```


I did not copy the code over here


## Run on efficiency

I probably need to change this because I swapped the order of the ENSG and gene in the gene annotation file above

```
coverageBed -sorted -counts -s -a Homo_sapiens.GRCh38.${release}.pause.bed -b ${name}_PE1_signal.bed -g hg38.chrom.order.txt | awk '$7>0' | sort -k5,5 -k7,7nr | sort -k5,5 -u > ${name}_pause.bed

#discard anything with chr and strand inconsistencies
join -1 5 -2 5 ${name}_pause.bed Homo_sapiens.GRCh38.${release}.bed | awk '{OFS="\t";} $2==$8 && $6==$12 {print $2, $3, $4, $1, $6, $7, $9, $10}' | awk '{OFS="\t";} $5 == "+" {print $1,$2+480,$8,$4,$6,$5} $5 == "-" {print $1,$7,$2 - 380,$4,$6,$5}' |  awk  '{OFS="\t";} $3>$2 {print $1,$2,$3,$4,$5,$6}' | sort -k1,1 -k2,2n > ${name}_pause_counts_body_coordinates.bed

#column ten is Pause index
coverageBed -sorted -counts -s -a ${name}_pause_counts_body_coordinates.bed -b ${name}_PE1_signal.bed -g hg38.chrom.order.txt | awk '$7>0' | awk '{OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,$5/100,$7/($3 - $2)}' | awk '{OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$8/$9}' > ${name}_pause_body.bed

pause_index.R ${name}_pause_body.bed
```

## Estimate nascent RNA purity with exon / intron density ratio
```
coverageBed -sorted -counts -s -a Homo_sapiens.GRCh38.${release}.introns.bed -b ${name}_PE1_signal.bed -g hg38.chrom.order.txt  | awk '$7>0' | awk '{OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,($3 - $2)}' > ${name}_intron_counts.bed

coverageBed -sorted -counts -s -a Homo_sapiens.GRCh38.${release}.no.first.exons.named.bed -b ${name}_PE1_signal.bed -g hg38.chrom.order.txt | awk '$7>0' | awk '{OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,($3 - $2)}' > ${name}_exon_counts.bed

exon_intron_ratio.R ${name}_exon_counts.bed ${name}_intron_counts.bed
```



