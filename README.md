# Nascent_RNA_Methods

# One time things
Get code (put them in a path or leave them in the current directory)
```
wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/insert_size.R
wget https://raw.githubusercontent.com/guertinlab/fqComplexity/main/complexity_pro.R
wget https://raw.githubusercontent.com/guertinlab/fqComplexity/main/fqComplexity
wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/pause_index.R
wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/exon_intron_ratio.R

chmod +x insert_size.R
chmod +x fqComplexity
chmod +x complexity_pro.R
chmod +x pause_index.R
chmod +x exon_intron_ratio.R
```

Build genomes
```
wget https://github.com/databio/ref_decoy/raw/master/human_rDNA.fa.gz
gunzip human_rDNA.fa.gz
bowtie2-build human_rDNA.fa human_rDNA

wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
bowtie2-build hg38.fa hg38
```

Parse gene annotations for use in calculating pause index and exon / intron density
```
#Go here to find the latest release number http://www.ensembl.org/index.html
release=104

wget http://ftp.ensembl.org/pub/release-${release}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${release}.chr.gtf.gz
gunzip Homo_sapiens.GRCh38.${release}.chr.gtf.gz

#don't know why the chr is not present
#parse all TSS--exons 1
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
    awk '{OFS="\t";} {print $1,$4,$5,$14,$10,$7}' | \
    sed 's/";//g' | \
    sed 's/"//g' | sed 's/chrMT/chrM/g' | sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.${release}.bed
    
sort -k4,4 Homo_sapiens.GRCh38.${release}.bed > Homo_sapiens.GRCh38.${release}.sorted.gene.bed

#identify and organize all exons within genes
intersectBed -s -a Homo_sapiens.GRCh38.${release}.bed -b Homo_sapiens.GRCh38.${release}.all.exons.bed | sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.${release}.all.exons.sorted.bed

#merge intervals that overlap
mergeBed -s -c 6 -o distinct -i Homo_sapiens.GRCh38.${release}.all.exons.sorted.bed | awk '{OFS="\t";} {print $1,$2,$3,$4,$2,$4}' | sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.${release}.all.exons.merged.bed

#remove all first exons (where pause is)
subtractBed -s -a Homo_sapiens.GRCh38.${release}.all.exons.merged.bed -b Homo_sapiens.GRCh38.${release}.tss.bed | sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.${release}.no.first.exons.bed


#all introns named
subtractBed -s -a Homo_sapiens.GRCh38.${release}.bed -b Homo_sapiens.GRCh38.${release}.all.exons.merged.bed | sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.${release}.introns.bed 

#get gene names of exons
intersectBed -s -wb -a Homo_sapiens.GRCh38.${release}.no.first.exons.bed -b Homo_sapiens.GRCh38.${release}.bed | awk '{OFS="\t";} {print $1,$2,$3,$10,$4,$4}' | sort -k1,1 -k2,2n >  Homo_sapiens.GRCh38.${release}.no.first.exons.named.bed

#then use coverageBed
#load into R and use aggregate for gene name

wget https://hgdownload-test.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
sort -k1,1 -k2,2n hg38.chrom.sizes | sed 's/chrMT/chrM/g' > hg38.chrom.order.txt
 
#first this
#window 20-120 
awk  '{OFS="\t";} $6 == "+" {print $1,$2+20,$2 + 120,$4,$5,$6} $6 == "-" {print $1,$3 - 120,$3 - 20,$4,$5,$6}' Homo_sapiens.GRCh38.${release}.tss.bed  | sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.${release}.pause.bed
```

# Initialize variables
```
release=104
read_size=30
UMI_length=8
cores=6
directory=/Users/guertinlab/Downloads/Batch1 

#For if you run this in a loop
name=$(echo $1 | awk -F"_PE1.fastq.gz" '{print $1}')

#For manual naming
name=LNCaP_10uMEnza_rep3_batch2
```

# Begin
```
cd $directory 
echo ${name}
gunzip ${name}_PE*.fastq.gz
```

First remove any bases beyond the adapter so that we can accurately identify duplicates
(I don't think we need to save the output text, since we get the insert size distribution later)
(Do you not want to do this in parallel?)
```
cutadapt --cores=0 --action=retain -a TGGAATTCTCGGGTGCCAAGG ${name}_PE1.fastq -o ${name}_PE1_trimmed.fastq
```

Deduplicate
```
fqdedup -i ${name}_PE1_trimmed.fastq -o ${name}_PE1_dedup.fastq
```

Remove adapters and filter adapter/adapter ligation products as well as 1 base inserts
I incorporate the fraction of reads that are adapter/adapter products as a metric below.
```
(cutadapt -m $((UMI_length+2)) -O 1 -a TGGAATTCTCGGGTGCCAAGG ${name}_PE1_dedup.fastq -o ${name}_PE1_noadap.fastq --too-short-output ${name}_PE1_short.fastq ) > ${name}_PE1_cutadapt.txt

(cutadapt -m $((UMI_length+2)) -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC ${name}_PE2.fastq -o ${name}_PE2_noadap.fastq --too-short-output ${name}_PE2_short.fastq ) > ${name}_PE2_cutadapt.txt

```

# DEGRADATION RNA INTEGRITY
```
reads=$(wc -l ${name}_PE1_noadap.fastq | awk '{print $1/4}')
fastq_pair -t $reads ${name}_PE1_noadap.fastq ${name}_PE2_noadap.fastq

flash -q --compress-prog=gzip --suffix=gz ${name}_PE1_noadap.fastq.paired.fq ${name}_PE2_noadap.fastq.paired.fq -o ${name}

./insert_size.R ${name}.hist ${UMI_length}

rm ${name}_PE*_noadap.fastq.paired.fq
```

# PROCESSING

Trim the UMI, remove inserts of less than 10 bases, and reverse complement
```
seqtk trimfq -b ${UMI_length} ${name}_PE1_noadap.fastq | seqtk seq -L 10 -r - > ${name}_PE1_processed.fastq
seqtk trimfq -e ${UMI_length} ${name}_PE2_noadap.fastq | seqtk seq -L 10 -r - > ${name}_PE2_processed.fastq
```

Remove reads aligning to rDNA
(I'm switching all flags to hex for consistency)
```
bowtie2 -p 3 -x human_rDNA -U ${name}_PE1_processed.fastq 2>${name}_bowtie2_rDNA.log | samtools sort -n - | samtools fastq -f 0x4 - > ${name}_PE1.rDNA.fastq

reads=$(wc -l ${name}_PE1.rDNA.fastq | awk '{print $1/4}')
fastq_pair -t $reads ${name}_PE1.rDNA.fastq ${name}_PE2_processed.fastq
```

Align to hg38
```
bowtie2 -p 3 -x hg38 --rf -1 ${name}_PE1.rDNA.fastq.paired.fq -2 ${name}_PE2_processed.fastq.paired.fq 2>${name}_bowtie2_hg38.log | samtools view -b - | samtools sort - -o ${name}.bam
```

# rDNA ALIGNMENT RATE
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

# FRACTION OF FILTERED READS THAT ARE MAPPABLE
this is a QC metric in itself. Concordant alignment rate is typically above 90% for good libraries
map is less stringently than concordantly mapped and most should map
```
map_pe1=$(samtools view -c -f 0x40 -F 0x4 ${name}.bam)
pre_alignment=$(wc -l ${name}_PE1.rDNA.fastq.paired.fq | awk '{print $1/4}')
alignment_rate=$(echo "scale=2 ; $map_pe1 / $pre_alignment" | bc)
```

# COMPLEXITY AND THEORETICAL READ DEPTH

calculate PE1 total raw reads
```
PE1_total=$(wc -l ${name}_PE1.fastq | awk '{print $1/4}')
```

calculate PE1 deduplicated reads
```
PE1_dedup=$(wc -l ${name}_PE1_dedup.fastq | awk '{print $1/4}')
```

calculate PE1 reads without adapters 
```
PE1_noadap=$(wc -l ${name}_PE1_noadap.fastq | awk '{print $1/4}')
```

This inverse of this factor is a QC metric for percent adapter/adapter ligation products (including 1 base inserts)

```
factorX=$(echo "scale=2 ; $PE1_dedup / $PE1_noadap" | bc)

echo fraction of reads that are not adapter/adapter ligation products
echo $factorX | awk '{print 1/$1}'
```

divide
```
factorY=$(echo "scale=2 ; $concordant_pe1 / $PE1_dedup" | bc)
```

two different curves:
this curve lets you know the quality of the library in terms of it's complexity.
if at 10 milion reads depth, 75% of the reads are unique, then it passes. The higher the better 

```
./fqComplexity -i ${name}_PE1_trimmed.fastq 
```

This curve is similar, but the goal is to estimate the raw read depth needed to acieve a target concordant aligned read count
the two factors needed are the fraction of the total reads that are adapter/adapter ligation products (this is the only step that preceeds deduplication)
 and the fraction of deduplicated reads that align concordantly after filtering rDNA-aligned reads, short reads, and unaligned reads
 simply solve for read_depth after providing a desired Concordant Aligned
rearrange the equation for them?
empirically noticed that copying and pasting the equation into R interprets the minus signs as hyphens

I removed factorX here based on deduplicating happening first (trimming the junk keeps all reads). Also this makes it run the whole thing again, which takes time. Should we instead just have a calculator that takes the output of the first run, factorY, and desired concordant aligned reads and outputs the number of raw reads needed?
```
./fqComplexity -i ${name}_PE1_trimmed.fastq -y $factorY
```

counterintuitively, you can have a high quality and complex library that is not practical to sequence to further depth because the number of adapter/adapter
reads is just too high
the post-depulication factors are all individual QC metrics, so these should be considered individually to determine whether the library is high quality.


#Get the reads in a BED
I changed the flags here
```
#samtools view -b -f 0x0040 ${name}.bam > ${name}_PE1.bam
samtools view -b -f 0x40 ${name}.bam > ${name}_PE1.bam
#samtools view -bh -F 20 ${name}_PE1.bam > ${name}_PE1_plus.bam
samtools view -bh -F 0x14 ${name}_PE1.bam > ${name}_PE1_plus.bam
samtools view -bh -f 0x10 ${name}_PE1.bam > ${name}_PE1_minus.bam
    
seqOutBias hg38.fa ${name}_PE1_plus.bam --no-scale --bed ${name}_PE1_plus.bed --bw=${name}_PE1_plus.bigWig --tail-edge --read-size=$read_size
seqOutBias hg38.fa ${name}_PE1_minus.bam --no-scale --bed ${name}_PE1_minus.bed --bw=${name}_PE1_minus.bigWig --tail-edge --read-size=$read_size

awk '{OFS="\t";} {print $1,$2,$3,$4,$5,"+"}' ${name}_PE1_plus.bed > ${name}_PE1_plus_strand.bed
awk '{OFS="\t";} {print $1,$2,$3,$4,$5,"-"}' ${name}_PE1_minus.bed > ${name}_PE1_minus_strand.bed

cat ${name}_PE1_plus_strand.bed ${name}_PE1_minus_strand.bed | sort -k1,1 -k2,2n > ${name}_PE1_signal.bed
```


I did not copy the code over here


# Run on efficiency

```
coverageBed -sorted -counts -s -a Homo_sapiens.GRCh38.${release}.pause.bed -b ${name}_PE1_signal.bed -g hg38.chrom.order.txt | awk '$7>0' | sort -k5,5 -k7,7nr | sort -k5,5 -u > ${name}_pause.bed

#discard anyhting with chr and strand inconsistencies
join -1 5 -2 4 ${name}_pause.bed Homo_sapiens.GRCh38.${release}.sorted.gene.bed | awk '{OFS="\t";} $2==$8 && $6==$12 {print $2, $3, $4, $1, $6, $7, $9, $10}' | awk '{OFS="\t";} $5 == "+" {print $1,$2+480,$8,$4,$6,$5} $5 == "-" {print $1,$7,$2 - 380,$4,$6,$5}' |  awk  '{OFS="\t";} $3>$2 {print $1,$2,$3,$4,$5,$6}' | sort -k1,1 -k2,2n > ${name}_pause_counts_body_coordinates.bed

#column ten is Pause index
coverageBed -sorted -counts -s -a ${name}_pause_counts_body_coordinates.bed -b ${name}_PE1_signal.bed -g hg38.chrom.order.txt | awk '$7>0' | awk '{OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,$5/100,$7/($3 - $2)}' | awk '{OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$8/$9}' > ${name}_pause_body.bed

./pause_index.R ${name}_pause_body.bed
```

# the next part is to estimate nascent RNA purity with exon / intron density ratio
```
coverageBed -sorted -counts -s -a Homo_sapiens.GRCh38.${release}.introns.bed -b ${name}_PE1_signal.bed -g hg38.chrom.order.txt  | awk '$7>0' | awk '{OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,($3 - $2)}' > ${name}_intron_counts.bed

coverageBed -sorted -counts -s -a Homo_sapiens.GRCh38.${release}.no.first.exons.named.bed -b ${name}_PE1_signal.bed -g hg38.chrom.order.txt | awk '$7>0' | awk '{OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,($3 - $2)}' > ${name}_exon_counts.bed

./exon_intron_ratio.R ${name}_exon_counts.bed ${name}_intron_counts.bed
```



