# Nascent_RNA_Methods

this repository will 

#this is just a placeholder
```
cd /Users/guertinlab/Downloads/Batch1 
```

#Could pass this in as argument 2 or just hard code it if it's unlikely to change
```
UMI_length=8

name=$(echo $1 | awk -F"_PE1.fastq.gz" '{print $1}')
name=LNCaP_10uMEnza_rep3_batch2
echo ${name}

gunzip ${name}_PE*.fastq.gz
```

```
(cutadapt -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG ${name}_PE1.fastq -o ${name}_PE1_noadap.fastq --too-short-output ${name}_PE1_short.fastq ) > ${name}_PE1_cutadapt.txt

(cutadapt -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC ${name}_PE2.fastq -o ${name}_PE2_noadap.fastq --too-short-output ${name}_PE2_short.fastq ) > ${name}_PE2_cutadapt.txt

```

```
fqdedup -i ${name}_PE1_noadap.fastq -o ${name}_PE1_dedup.fastq

```

# DEGRADATION RNA INTEGRITY

```
reads=$(wc -l ${name}_PE1_dedup.fastq | awk '{print $1/4}')
fastq_pair -t $reads ${name}_PE1_dedup.fastq ${name}_PE2_noadap.fastq

flash -q --compress-prog=gzip --suffix=gz ${name}_PE1_dedup.fastq.paired.fq ${name}_PE2_noadap.fastq.paired.fq -o ${name}

Rscript ~/rscripts/insert_size.R ${name}.hist ${UMI_length}

#I think we can delete these they are only used for insert size 
rm ${name}_PE2_noadap.fastq.paired.fq
rm ${name}_PE1_dedup.fastq.paired.fq

```

# PROCESSING

```
seqtk trimfq -b ${UMI_length} ${name}_PE1_dedup.fastq | seqtk seq -L 10 -r - > ${name}_PE1_processed.fastq

wget https://github.com/databio/ref_decoy/raw/master/human_rDNA.fa.gz

gunzip human_rDNA.fa.gz
bowtie2-build human_rDNA.fa human_rDNA
bowtie2 -p 3 -x human_rDNA -U ${name}_PE1_processed.fastq 2>${name}_bowtie2_rDNA.log | samtools sort -n - | samtools fastq -f 4 - > ${name}_PE1.rDNA.fastq

seqtk trimfq -e ${UMI_length} ${name}_PE2_noadap.fastq | seqtk seq -L 10 -r - > ${name}_PE2_noadap_trimmed.fastq


reads=$(wc -l ${name}_PE1.rDNA.fastq | awk '{print $1/4}')
fastq_pair -t $reads ${name}_PE1.rDNA.fastq ${name}_PE2_noadap_trimmed.fastq
```

```
#
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
bowtie2-build hg38.fa hg38

bowtie2 -p 3 -x hg38 --rf -1 ${name}_PE1.rDNA.fastq.paired.fq -2 ${name}_PE2_noadap_trimmed.fastq.paired.fq 2>${name}_bowtie2_hg38.log | samtools view -b - | samtools sort - -o ${name}.bam
```




# rDNA ALIGNMENT RATE

```
PE1_prior_rDNA=$(wc -l ${name}_PE1_processed.fastq | awk '{print $1/4}')
PE1_post_rDNA=$(wc -l ${name}_PE1.rDNA.fastq | awk '{print $1/4}')
```

rRNA alignment rate (does not account for low genome alignment rates, so this can be artifically low if the concordant alignment rates are low)
if concodarnt alignmnet rate are low this supercedes rDNA alignment rate and multiplying by the inverse of the concordant alugnment rate gives a better approximation of rDNA alignmetn rate
```
not_considering_overall_alignment_rate=$(echo "$(($PE1_prior_rDNA-$PE1_post_rDNA))" | awk -v myvar=$PE1_prior_rDNA '{print $1/myvar}')
```
alternatively, of the aligned reads, what fractino is rDNA:
this is what PEPRO should do
```
overall_alignment_considered=$(echo "$(($PE1_prior_rDNA-$PE1_post_rDNA))" | awk -v myvar=$concordant_pe1 '{print $1/myvar}')
```




# FRACTION OF FILTERED READS THAT ARE MAPPABLE
this is a QC metric in itself. Concordant alignment rate is typically above 90% for good libraries
map is less stringently than concordantly mapped and most should map

```
map_pe1=$(samtools view -c -f 0x40 -F 0x4 ${name}.bam)
pre_alignment=$(wc -l ${name}_PE1.rDNA.fastq.paired.fq | awk '{print $1/4}')
alignment_rate=$(echo "scale=2 ; $map_pe1 / $pre_alignment" | bc)
```



# COMPLEXITY AND THEORETICAL READ DEPTH

```
wget https://raw.githubusercontent.com/guertinlab/fqComplexity/main/complexity_pro.R
wget https://raw.githubusercontent.com/guertinlab/fqComplexity/main/fqComplexity
chmod +x fqComplexity
chmod +x complexity_pro.R
```
put them in a path or leave them in the current directory


calculate PE1 total raw reads
```
PE1_total=$(wc -l ${name}_PE1.fastq | awk '{print $1/4}')
```

calculate PE1 reads without adapters 
```
PE1_noadap=$(wc -l ${name}_PE1_noadap.fastq | awk '{print $1/4}')
```
divide
```
factorX=$(echo "scale=2 ; $PE1_total / $PE1_noadap" | bc)
```
calculate dedup PE1 reads 
```
PE1_dedup=$(wc -l ${name}_PE1_dedup.fastq | awk '{print $1/4}')
```
need to extract aligned reds from ${name}_bowtie2_hg38.log
```
concordant_pe1=$(samtools view -c -f 0x42 ${name}.bam)
```
divide
```
factorY=$(echo "scale=2 ; $concordant_pe1 / $PE1_dedup" | bc)
```

two different curves:
this curve lets you know the quality of the library in terms of it's complexity.
if at 10 milion reads depth, 75% of the reads are unique, then it passes. The higher the better 

```
./fqComplexity -i ${name}_PE1_noadap.fastq 
```

This curve is similar, but the goal is to estimate the raw read depth needed to acieve a target concordant aligned read count
the two factors needed are the fraction of the total reads that are adapter/adapter ligation products (this is the only step that preceeds deduplication)
 and the fraction of deduplicated reads that align concordantly after filtering rDNA-aligned reads, short reads, and unaligned reads
 simply solve for read_depth after providing a desired Concordant Aligned
rearrange the equation for them?
empirically noticed that copying and pasting the equation into R interprets the minus signs as hyphens

```
./fqComplexity -i ${name}_PE1_noadap.fastq -x $factorX -y $factorY
```

counterintuitively, you can have a high quality and complex library that is not practical to sequence to further depth because the number of adapter/adapter
reads is just too high
the post-depulication factors are all individual QC metrics, so these should be considered individually to determine whether the library is high quality.


# get the reads in a BED

```
samtools view -b -f 0x0040 ${name}.bam > ${name}_PE1.bam
samtools view -bh -F 20 ${name}_PE1.bam > ${name}_PE1_plus.bam
samtools view -bh -f 0x10 ${name}_PE1.bam > ${name}_PE1_minus.bam
    
seqOutBias hg38.fa ${name}_PE1_plus.bam --no-scale --bed ${name}_PE1_plus.bed --bw=${name}_PE1_plus.bigWig --tail-edge --read-size=30
seqOutBias hg38.fa ${name}_PE1_minus.bam --no-scale --bed ${name}_PE1_minus.bed --bw=${name}_PE1_minus.bigWig --tail-edge --read-size=30
```

# the next part is to estimate nascent RNA purity with exon / intron density ratio
```
wget http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gtf.gz
#wget ftp://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh38.87.gtf.gz
gunzip Homo_sapiens.GRCh38.104.chr.gtf.gz

#parse all TSS--exons 1
grep 'exon_number "1"' Homo_sapiens.GRCh38.104.chr.gtf | \
    sed 's/^/chr/' | \
    awk '{OFS="\t";} {print $1,$4,$5,$14,$20,$7}' | \
    sed 's/";//g' | \
    sed 's/"//g' | sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.104.tss.bed

#extract all exons
grep 'exon_number' Homo_sapiens.GRCh38.104.chr.gtf | \
    sed 's/^/chr/' | \
    awk '{OFS="\t";} {print $1,$4,$5,$14,$20,$7}' | \
    sed 's/";//g' | \
    sed 's/"//g' | sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.104.all.exons.bed


#extract all complete gene annotations
awk '$3 == "gene"' Homo_sapiens.GRCh38.104.chr.gtf | \
    sed 's/^/chr/' | \
    awk '{OFS="\t";} {print $1,$4,$5,$14,$10,$7}' | \
    sed 's/";//g' | \
    sed 's/"//g' | sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.104.bed


#identify and organize all exons within genes
intersectBed -s -a Homo_sapiens.GRCh38.104.bed -b Homo_sapiens.GRCh38.104.all.exons.bed | sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.104.all.exons.sorted.bed

#merge intervals that overlap
mergeBed -s -c 6 -o distinct -i Homo_sapiens.GRCh38.104.all.exons.sorted.bed | awk '{OFS="\t";} {print $1,$2,$3,$4,$2,$4}' | sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.104.all.exons.merged.bed

#remove all first exons (where pause is)
subtractBed -s -a Homo_sapiens.GRCh38.104.all.exons.merged.bed -b Homo_sapiens.GRCh38.104.tss.bed | sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.104.no.first.exons.bed


#all introns named
subtractBed -s -a Homo_sapiens.GRCh38.104.bed -b Homo_sapiens.GRCh38.104.all.exons.merged.bed | sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.104.introns.bed 

#get gene names of exons
intersectBed -s -wb -a Homo_sapiens.GRCh38.104.no.first.exons.bed -b Homo_sapiens.GRCh38.104.bed | awk '{OFS="\t";} {print $1,$2,$3,$10,$4,$4}' | sort -k1,1 -k2,2n >  Homo_sapiens.GRCh38.104.no.first.exons.named.bed


```

