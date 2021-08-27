# Nascent_RNA_Methods

Precision genomic run-on assays (PRO-seq) quantify nascent RNA at single nucleotide resolution with strand specificity. Here, we deconstruct a recently published genomic nascent RNA processing pipeline (PEPPRO) into its components. The analyses are presented as individual code chunks with comprehensive details. PRO-seq experiments are evolving and variations can be found throughout the literature. We describe analysis of paired end libraries that contain a unique molecular identifier on the the 5´ adapter adjacent to the RNA insert. These two features allow for calculations of RNA fragment length and library complexity. Our hope is that the end users will adopt the analysis framework and modify the workflow as needed. We describe and provide the code to calculate the following quality control metrics: library complexity, nascent RNA purity, nuclear run on efficiency, alignment rate, sequencing depth, and RNA integrity.        

# One time software installations, downloads, and processing steps

Many processes, downloads, and software installations are reused in analyses. Although we refer to these as "one-time" steps, annotations are often updated and newer software versions are released.  

## Dependencies, software, and scripts

We present novel quality control metrics and specialized software herein, but much of the workflow depends upon more general software. Fortunately, this software is well-maintained and documented, so we only provide a short description and the links below. 

[Cite all the software in line]

`bedtools` :  a comprehensive suite of tools that efficiently perform a wide range of operations on genomic intervals.  https://bedtools.readthedocs.io/en/latest/

`bowtie2` : aligns sequence reads to reference sequences. http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

`cutadapt` : removes a defined sequence, such as adapter sequence, from sequence reads.  https://cutadapt.readthedocs.io/en/stable/

`fastq_pair` : outputs only sequence reads that have a matched paired end read and a separates unpaired reads. https://github.com/linsalrob/fastq-pair 

`FLASH` : merges paired end reads by detecting overlapping sequence. https://ccb.jhu.edu/software/FLASH/ 

`fqdedup` : removes duplicated sequences from FASTQ files with a small memory footprint. https://github.com/guertinlab/fqdedup

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

The goal of the following operations is to output a set of exons that excludes all instances of first exons, output all introns, and output a set of all potential pause regions for a gene by taking the region for 20 - 120 downstream of all exon 1 annotations. The `mergeBed` command collapses all overlapping intervals and the gene name information is lost. We exclude all first exon coordinates with `subtractBed`, then `intersectBed` reassigns gene names to all remaining exons. The final `awk` command defines a 100 base pause region window downstream of all transcription start sites based on the gene strand.    

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

PRO-seq data can be analyzed in many sophisticated ways, including defining primary transcripts (cite primaryTranscriptAnnotation), identifying putative enhancers (cite dREG), detecting prominent transcription start sites (cite something I saw come out recently but haven't tested yet), or quantifying changes in transcription (cite Fabiana and Leighton). Here, we only describe quality control metrics that are used to determine if a PRO-seq library is worth analyzing in depth. 

## Initialize variables

The naming convention we recommend is the following: cellType_conditions_replicate_pairedend.fastq.gz. For example, a gzipped  paired end 1 (PE1) file from the second replicate of treating MCF7 cells with estrogen (E2) for 20 minutes would be: `MCF7_20minE2_rep2_PE1.fastq.gz`. Many of the lines of code assume this naming convention, especially with regards to the trailing `_PE1.fastq.gz` and `_PE2.fastq.gz`. 

We initialize six variables at the start:

`$directory`: location of the sequencing file

`$filename`: name of the gzipped FASTQ file

`$release`: Ensembl gene annotation release number

`$UMI_length`: length of the UMI on the 5´ end of the paired end 1 read.

`$read_size`: read length minus UMI length.

`$cores`: number of cores for parallel processing

```
directory=/Users/guertinlab/sequencing_run1 
filename=MCF7_20minE2_rep2_PE1.fastq.gz
release=104
UMI_length=8
read_size=30
cores=6
```

## Preprocessing 

Navigate to the folder that contains the sequencing file and gunzip the file.

```
cd $directory 
name=$(echo $filename | awk -F"_PE1.fastq.gz" '{print $1}')
echo $name
gunzip ${name}_PE*.fastq.gz
```

## Processing reads 

The first processing step is to remove adapter sequence and we simultaneously discard reads that have insert sizes of one base. The 3´ adapter contains a UMI, which is sequenced prior to the adapter. Therefore, the vast majority of adapter/adapter ligation products have read lengths of exactly the length of the UMI. The option `-m $((UMI_length+2))` provides a one base buffer and discards reads with a length of the UMI + 1.

The fraction of reads that result from adapter/adapter ligation products can be a useful metric. We calculate this value by first counting the number of lines in the original FASTQ file using `wc -l` and divide that value by 4 using `awk '{print $1/4}'` because FASTQ files contain four lines per sequence entry. This same operation is performed on the output file of reads that were too short, in this case 0 or 1 base insert. Divide the adapter/adapter ligation product value by the total and round to the hundredth with `$(echo "scale=2 ; $PE1_w_Adapter / $PE1_total" | bc)`.

This value varies widely depending upon whether a size selection was performed in the library preparation. We recently dropped the size selection step from the protocol in an effort to reduce bias against small RNA inserts (cite Sathyan). If no size selection is performed, this value can be quite high; however, a high value does not indicate poor library quality. Because sequencing is relatively cheap, we tolerate up to 80% adapter/adapter ligation products. One needs to balance to cost of performing another experiment with sequencing uninformative adapter sequences. Later on we provide a formula for determining the required sequncing depth to result in a desired number of concordant aligned reads. So why this is a useful number, if all other QC metrics point to high quality data, then we recommend further sequencing depth if less than 80% of the reads are adapter/adapter ligation products. 

Each quality control metric can be distilled down to a single number that is printed to the file `${name}_QC_metrics.txt`. Alongside the value, we include the recommended threshold. We will continue to print all metrics to this file and plot the data at the end of the workflow.


```
cutadapt --cores=$cores -m $((UMI_length+2)) -O 1 -a TGGAATTCTCGGGTGCCAAGG ${name}_PE1.fastq -o ${name}_PE1_noadap.fastq --too-short-output ${name}_PE1_short.fastq > ${name}_PE1_cutadapt.txt
cutadapt --cores=$cores -m $((UMI_length+10)) -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC ${name}_PE2.fastq -o ${name}_PE2_noadap.fastq --too-short-output ${name}_PE2_short.fastq > ${name}_PE2_cutadapt.txt

PE1_total=$(wc -l ${name}_PE1.fastq | awk '{print $1/4}')
PE1_w_Adapter=$(wc -l ${name}_PE1_short.fastq | awk '{print $1/4}')
AAligation=$(echo "scale=2 ; $PE1_w_Adapter / $PE1_total" | bc)

echo -e  "value\texperiment\tthreshold\tmetric" > ${name}_QC_metrics.txt
echo -e "$AAligation\t$name\t0.80\tAdapter/Adapter" >> ${name}_QC_metrics.txt
```

The next step removes reads with a length shorter than 10 bases and reverse complements the file so that the aligned read corresponds to the appropriate reference genome strand.

```
seqtk seq -L $((UMI_length+10)) -r ${name}_PE1_noadap.fastq > ${name}_PE1_noadap_trimmed.fastq 
```

PRO-seq can have several independent reads that have the same genomic ends because promoter proximal pausing positions can be focused and the 5´end of the RNA is often the transcription start site. One cannot filter potential PCR duplicate reads based on whether two independent reads are identical as determined by having identical paired end reads alignments. Therefore, we remove duplicate sequences from the FASTQ file based on the presence of the UMI. By pairing the PE1 and PE2 reads, we effectively deduplicate the PE2 based on the presence of the PE1 UMI. 

```
fqdedup -i ${name}_PE1_noadap_trimmed.fastq -o ${name}_PE1_dedup.fastq

#this variable is a near-optimal optimal table size value for fastq_pair: 
PE1_noAdapter=$(wc -l ${name}_PE1_noadap.fastq | awk '{print $1/4}')

fastq_pair -t $PE1_noAdapter ${name}_PE1_noadap.fastq ${name}_PE2_noadap.fastq
```

## RNA integrity score

We measure RNA degradation by searching for overlap between paired end reads with `flash` and the resultant histogram output. We empirically found that there are fewer reads within the range of 10 - 20 than the range of 30 - 40 for high quality libraries. RNA only starts to protrude from the RNA Polymerase II exit channel at approximately 20 bases in length, so 20 base of RNA is protected from degradation during the run on. Libraries with a substantial amount of degradation after the run on step are enriched for species in the range 10 - 20. A degradation ratio of less than 1 indicates a high quality library.  

```
flash -q --compress-prog=gzip --suffix=gz ${name}_PE1_noadap.fastq.paired.fq ${name}_PE2_noadap.fastq.paired.fq -o ${name}
insert_size.R ${name}.hist ${UMI_length}

rm ${name}_PE*_noadap.fastq.paired.fq
```

## Processing for alignment

The final processing step removes the UMI from both paired end reads and reverse complements the paired end 2 read.

```
seqtk trimfq -e ${UMI_length} ${name}_PE1_noadap_trimmed.fastq  > ${name}_PE1_processed.fastq
seqtk trimfq -e ${UMI_length} ${name}_PE2_noadap.fastq | seqtk seq -r - > ${name}_PE2_processed.fastq
```

## Remove reads aligning to rDNA

By first aligning to the rDNA we can later estimate nascent RNA purity and avoid spurious read pile ups at region in the genome that non-uniquely align to both the rDNA locus and elsewhere in the genome. While between 70 - 80% of stable RNA is rRNA, generally between 10 - 15% of the nascent RNA arises from rRNA. Even 10% of the library aligning to rDNA loci is extremely enriched, so any reads that map non-uniquely to both rDNA and non-rDNA regions in the genome result in huge artifactual spikes in the data if rDNA-aligned reads are not first removed. The `   
By first aligning to the rDNA we can later estimate nascent RNA purity and avoid spurious read pile ups at region in the genome that non-uniquely align to both the rDNA locus and elsewhere in the genome. While between 70 - 80% of stable RNA is rRNA, generally between 10 - 15% of the nascent RNA arises from rRNA. Even 10% of the library aligning to rDNA loci is extremely enriched, so any reads that map non-uniquely to both rDNA and non-rDNA regions in the genome result in huge artifactual spikes in the data if rDNA-aligned reads are not first removed. The `-f 0x4` flag of the `samtools fastq` command specifies that only unmapped reads should be included in the FASTQ output.    

```
bowtie2 -p $cores --maxins 1000 -x human_rDNA -U ${name}_PE1_processed.fastq 2>${name}_bowtie2_rDNA.log | samtools sort -n - | samtools fastq -f 0x4 - > ${name}_PE1.rDNA.fastq

#this effectively removes PE2-aligned reads with a rDNA-aligned mate
reads=$(wc -l ${name}_PE1.rDNA.fastq | awk '{print $1/4}')
fastq_pair -t $reads ${name}_PE1.rDNA.fastq ${name}_PE2_processed.fastq
```

## Genome alignment
The last processing step for individual libraries is to align to genome. Note that we reverse complemented both reads, so the `--rf` flag indicates that the reads are oriented opposite to typical RNA-seq data. The `samtools` commands convert the file to a compressed binary BAM format and sort the reads. 

```
bowtie2 -p cores --maxins 1000 -x hg38 --rf -1 ${name}_PE1.rDNA.fastq.paired.fq -2 ${name}_PE2_processed.fastq.paired.fq 2>${name}_bowtie2_hg38.log | samtools view -b - | samtools sort - -o ${name}.bam
```

## rDNA alignment rate
<!--- ```
%% PE1_prior_rDNA=$(wc -l ${name}_PE1_processed.fastq | awk '{print $1/4}')
%% PE1_post_rDNA=$(wc -l ${name}_PE1.rDNA.fastq | awk '{print $1/4}')
%% ```


rRNA alignment rate (does not account for low genome alignment rates, so this can be artifically low if the concordant alignment rates are low)
If concodarnt alignmnet rate are low this supercedes rDNA alignment rate and multiplying by the inverse of the concordant alignment rate gives a better approximation of rDNA alignment rate

```
not_considering_overall_alignment_rate=$(echo "$(($PE1_prior_rDNA-$PE1_post_rDNA))" | awk -v myvar=$PE1_prior_rDNA '{print $1/myvar}')
```
--->
In order to rDNA alignment rate, we first counts the number of reads prior to rDNA alignment and after removing rDNA aligned reads. By subtracting the post-alignment read count from the input, we calculate the total number of rDNA-aligned reads. The command `samtools view -c -f 0x42` counts the concordantly aligned paired end 1 reads. Finally we calcualte the fraction of aligned reads that map to the rDNA genome and print it to the QC metrics file.    

```
PE1_prior_rDNA=$(wc -l ${name}_PE1_processed.fastq | awk '{print $1/4}')
PE1_post_rDNA=$(wc -l ${name}_PE1.rDNA.fastq | awk '{print $1/4}')
total_rDNA=$(echo "$(($PE1_prior_rDNA-$PE1_post_rDNA))") 

concordant_pe1=$(samtools view -c -f 0x42 ${name}.bam)
total_concordant=$(echo "$(($concordant_pe1+$total_rDNA))")

rDNA_alignment=$(echo "scale=2 ; $total_rDNA / $total_concordant" | bc)

echo -e "$rDNA_alignment\t$name\t0.20\trDNA Alignment Rate" >> ${name}_QC_metrics.txt

```

## Mappability rate
After all the processing, the vast majority of reads should map concordantly to the genome. Concordant alignment rate for successful PRO-seq experiments is typically above 90%. Again, we count specific reads using `samtools` and by counts reads in FASTQ files, then calcuate alignment rate. We recommend the followign site to help understand the meaning of samtools flags: https://broadinstitute.github.io/picard/explain-flags.html.

```
map_pe1=$(samtools view -c -f 0x40 -F 0x4 ${name}.bam)
pre_alignment=$(wc -l ${name}_PE1.rDNA.fastq.paired.fq | awk '{print $1/4}')
alignment_rate=$(echo "scale=2 ; $map_pe1 / $pre_alignment" | bc)

echo -e "$alignment_rate\t$name\t0.90\tAlignment Rate" >> ${name}_QC_metrics.txt
```

## Complexity and theoretical read depth

We developed `fqComplexity` to serve two purposes: 1) calculate the number of reads that are non-PCR duplicates as a metric for complexity; and 2) provide a formula and constants to calculate the theoretical read depth that will result in a user-defined number of concordant aligned reads. The equation accounts for all upstream processing steps. Over 10 milion concordantly aligned reads is typically sufficient if you have 3 or more replicates. 

For the first metrics, the input FASTQ file is preprocessed and the UMI is still included in the FASTQ DNA sequence. The FASTQ file is subsampled into deciles and the intermediate files are deduplicated. The input and output numbers are logged in `${name}_complexity.log`. An asymptotic regression model is fit to the data and the total number of unique reads at 10 million read depth is printed on the resulting PDF plot. We recommend that at least 7.5 million reads are unique at a depth of 10 million (i.e. 75% of reads are unique if you were to align exactly 10 million reads).

```
fqComplexity -i ${name}_PE1_noadap_trimmed.fastq
```

We need two factors to derive the constants to calculate the theoretical read depth for a specified concordant aligned read depth. First, we divide the PE1 total raw reads by processed PE1 reads that have adapter contamination and that have inserts less than 10 bases removed. The second value is the fraction of deduplicated reads that align concordantly to the non-rDNA genome.  

```
PE1_total=$(wc -l ${name}_PE1.fastq | awk '{print $1/4}')
PE1_noadap_trimmed=$(wc -l ${name}_PE1_noadap_trimmed.fastq | awk '{print $1/4}')

factorX=$(echo "scale=2 ; $PE1_total / $PE1_noadap_trimmed" | bc)

echo fraction of reads that are not adapter/adapter ligation products or below 10 base inserts
echo $factorX | awk '{print 1/$1}'

#calculate PE1 deduplicated reads
PE1_dedup=$(wc -l ${name}_PE1_dedup.fastq | awk '{print $1/4}')

#divide
factorY=$(echo "scale=2 ; $concordant_pe1 / $PE1_dedup" | bc)
```

Invoke `fqComplexity` with the `-x` and `-y` options specified to fit an asymptotic regression model to the factor-scaled log file.

```
fqComplexity -i ${name}_PE1_noadap_trimmed.fastq -x $factorX -y $factorY
```

Counterintuitively, you can have a high quality and complex library that is not practical to sequence to further depth because the number of adapter/adapter
reads is too damn high. The QC metrics should be considered to determine whether the library is high quality, but if the library is deemed high quality but you have low sequencing depth use this equation.  


## Get the reads in a BED
First we extract all the paired end 1 reads and separate reads based on their alignment to the plus or minus strand.

The software `seqOutBias` was originally developed to correct sequence bias from molecular genomics data. Although we are not correcting enzymatic sequence bias herein, there are many features of `seqOutBias` that are useful. Note that we include the `--no-scale` option to not correct sequence bias. The software outputs a bigWig and BED file, but it also calculates mappability at the specified read length and excludes non-uniquely mappable reads. Lastly, invoking `--tail-edge` realigned the end of the read so that the exact position of RNA Polymerase is specified in the output BED and bigWig files.

The `awk` command turns the respective stranded files into BED6 files with strand information in column 6, then the files are concatenated and sorted.
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


## Run on efficiency

RNA polymerases that are associated with gene bodies efficiently incorporate numcleotides during the run on reaction under most conditions, but promoter proximal paused RNA polymerase require high salt or detergent to run on efficiently (cite Leighton cells reports and ROugvie 1988). Therefore, the pause index is used to quantify run on efficiency. Pause index is the density of signal in the promoter-proximal pause region divided by density in the gene body. Hoever, since pause windows are variable, pause indices can differ substantially depending upon what one considers the pause window. There are many exon 1 gene annotations depending on gene isoforms and the upstream most annotated TSS is not necessarily the prominently transcribed isoform. It is common practice to choose the upstream most TSS, but this will cause the pause index to be deflated. Here, we define the pause window for a gene as position 20 -120 downstream of the most prominent TSS. The most prominent TSS is determined by calcualating the density in this 20 - 120 window for all annotated TSSs and choosing the TSS upstream of the most RNA-polymerase dense region.   

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



