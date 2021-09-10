---
title: "Processing and evaluating the quality of genome-wide nascent transcription profiling libraries"
author:
- "Thomas G. Scott"
- "André L. Martins"
- "Michael J. Guertin"

---


# Abstract

Precision genomic run-on assays (PRO-seq) quantify nascent RNA at single nucleotide resolution with strand specificity. Here we deconstruct a recently published genomic nascent RNA processing pipeline (PEPPRO) into its components and link the analyses to the underlying molecular biology. PRO-seq experiments are evolving and variations can be found throughout the literature. The analyses are presented as individual code chunks with comprehensive details so that users can modify the framework to accommodate different protocols. We present the framework to quantify the following quality control metrics: library complexity, nascent RNA purity, nuclear run on efficiency, alignment rate, sequencing depth, and RNA integrity.        


# Nature of PRO-seq experiments

Genomic nascent RNA profiling, such as precision genomic run-on assays (PRO-seq) [@kwak2013precise], quantify the precise position and direction of transcriptionally engaged RNA polymerases in the genome. Nascent RNA profiling methods have many advantages compared to measuring RNA levels with conventional RNA-seq. For instance, PRO-seq directly measures active transcription, compared to RNA-seq which quantifies steady-state RNA levels that are influenced by RNA stability. In fact, one can leverage the discordance between RNA-seq and PRO-seq expression to estimate genome-wide RNA half-lives [@blumberg2021characterizing]. Unstable bidirectional transcripts are a hallmark of enhancers and promoters and these unstable transcripts are used to directly infer regulatory element activity [@wang2019identification; @scruggs2015bidirectional]. Additionally, because PRO-seq measures active transcription it is used to detect immediate changes in transcription in time-course data need for RNAs to accumulate or degrade. Lasly, nascent RNA profiling can determines RNA polymerase location in the genome and within any genomic feature. Changes in RNA polymerase distribution can inform on how various treatments and stimuli regulate steps in the transcription cycle [@sathyan2019improved; @danko2013signaling; @hah2011rapid]. New genomic nacent RNA-seq methodologies [@schwalb2016tt; @scruggs2015bidirectional; @churchman2011nascent; @muhar2018slam; @mayer2015native] and increasing nascent RNA profiling data production necessitate flexible analysis workflows and standardized quality control metrics [@smith2021peppro]. Here we describe processing and analysis of paired end PRO-seq libraries with unique molecular identifiers presented as deconstructed code that can be easily adapted to fit a diversity of protocols and experimental details. 


# One time software installations, downloads, and processing steps

Many processes, downloads, and software installations are reused throughout the analyses. Although we refer to these as "one-time" steps, users should periodically check for updated annotations and new software releases. 

## Dependencies, software, and scripts

We present specialized software and scripts herein, but much of the workflow depends upon more general software. These general bioinformatic software tools are well-maintained and documented, so we provide short descriptions and the links below. 

`bedtools` :  a comprehensive suite of tools that efficiently perform a wide range of operations on genomic intervals.  https://bedtools.readthedocs.io/en/latest/ [@quinlan2010bedtools]

`bowtie2` : aligns sequencing reads to reference sequences. http://bowtie-bio.sourceforge.net/bowtie2/index.shtml [@langmead2012fast]

`cutadapt` : removes a defined sequence, such as adapter sequence, from sequencing reads.  https://cutadapt.readthedocs.io/en/stable/ [@martin2011cutadapt]

`fastq_pair` : outputs only sequencing reads that have a matched paired end read. https://github.com/linsalrob/fastq-pair [@edwards2019fastq]

`FLASH` : merges paired end reads by detecting overlapping sequence. https://ccb.jhu.edu/software/FLASH/ [@magovc2011flash] 

`fqdedup` : removes duplicated sequences from FASTQ files. https://github.com/guertinlab/fqdedup [@martins2018fqdedup]

`samtools` : a suite of tools for parsing and interfacing with high throughput sequencing data files. http://www.htslib.org [@li2009sequence]

`seqOutBias` : software that parses files and outputs desired formats with the option to correct enzymatic sequence biases. https://github.com/guertinlab/seqOutBias [@martins2018universal]

`seqtk` : a multifunctional toolkit for processing sequence files, including trimming a defined number of bases from the ends of reads and reverse complementing sequencing reads. https://github.com/lh3/seqtk [@li2013seqtk]

`wget` : retrieves files from a wide range of internet protocols. https://www.gnu.org/software/wget/

`R` packages:

`lattice`: graphics plotting package.  https://cran.r-project.org/web/packages/lattice/lattice.pdf [@sarkar2008lattice]

`DESeq2`: statistical package for quantifying differences in counts-based genomics data. https://bioconductor.org/packages/release/bioc/html/DESeq2.html [@love2014moderated]

In addition, we developed the following software and R scripts to facilitate data analysis and graphical output. Below, we use `wget` to retrieve the software and scripts. The command `chmod +x` changes the permissions of the files to executable. 

\scriptsize
```bash 
wget https://raw.githubusercontent.com/guertinlab/fqComplexity/main/fqComplexity
wget https://raw.githubusercontent.com/guertinlab/fqComplexity/main/complexity_pro.R
wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/insert_size.R
wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/pause_index.R
wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/exon_intron_ratio.R
wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/plot_all_metrics.R
wget https://raw.githubusercontent.com/guertinlab/Nascent_RNA_Methods/main/differential_expression.R

chmod +x insert_size.R
chmod +x fqComplexity
chmod +x complexity_pro.R
chmod +x pause_index.R
chmod +x exon_intron_ratio.R
chmod +x plot_all_metrics.R
chmod +x differential_expression.R
```
\normalsize

Next, move the software dependencies and R scripts to a directory within the `$PATH` variable.      

## Reference genomes

PRO-seq experiments have been performed in a variety of organisms including yeast [@booth2016divergence], Drosophila [@kwak2013precise; @duarte2016transcription], and humans [@core2014analysis]. Analysis of the data requires alignment to a reference genome annotation. The first step is to use `wget` to retrieve the reference genome. Many websites host the assembly data in FASTA format, such as the human genome build 38 shown below retrieved from the UCSC genome browser server [@karolchik2003ucsc]. The `gunzip` command unzips the reference genome file and `bowtie2-build` indexes the file to allow for efficient alignment. The code also retrieves, unzips, and builds the human rDNA reference genome [@stolarczyk2020refgenie] so that we can calculate rDNA alignment rates as a metric for nascent RNA purity.

```
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
bowtie2-build hg38.fa hg38

wget https://github.com/databio/ref_decoy/raw/master/human_rDNA.fa.gz
gunzip human_rDNA.fa.gz
bowtie2-build human_rDNA.fa human_rDNA
```

## Reference gene annotation

The quality control metrics outlined herein require the counting of sequence reads that align to three genomic features: exons, intron, and promoter-proximal pause regions. Gene annotations are available from many sources and we outline retrieval and parsing of GTF files from Ensembl [@yates2020ensembl]. The Ensembl website (http://www.ensembl.org/index.html) contains the information for the latest release, which at the time of writing this manuscript is release 104 for hg38. After retrieving and unzipping the file we parse out all exon 1 annotations. These coordinates include all annotated transcription start sites that correspond to different gene isoforms. Ensembl chromosome numbers do not include the preceding "chr", so the first `sed` command appends "chr" to the chromosome name for downstream compatibility with the reference genome chromosome names. The output is then piped to `awk`, which prints the following fields: chromosome coordinates in columns 1-3, gene name, Ensembl transcript ID (ENST), and strand.  Subsequent `sed` commands drop the semicolon and quote characters from the gene and Ensembl IDs while editing the mitochondrial chromosome to match the reference genome, "chrM" replacing "chrMT". Finally the exon output is sorted by the first column (chromosome), then the second column (starting position), in ascending order. The gene annotations are processed similarly, except the Ensembl gene ID (ENSG) replaces ENST in column 5.        

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
    sed 's/"//g' | sed 's/chrMT/chrM/g' | sort -k1,1 -k2,2n > Homo_sapiens.GRCh38.${release}.bed
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

PRO-seq data can be analyzed in many sophisticated ways, including defining primary transcripts [anderson2020defining], identifying putative enhancers [@wang2019identification], detecting prominent transcription start sites [@zhao2021deconvolution], or quantifying changes in RNA Polymerase density in different genic features [@core2012defining; duarte2016transcription; sathyan2019improved]. Here, we only describe quality control metrics that are used to determine if a PRO-seq library is worth analyzing in depth. 

## Initialize variables

The naming convention we recommend is the following: cellType_conditions_replicate_pairedend.fastq.gz. For example, a gzipped  paired end 1 (PE1) file from the second replicate of treating MCF7 cells with estrogen (E2) for 20 minutes would be: `MCF7_20minE2_rep2_PE1.fastq.gz`. Many of the lines of code assume this naming convention, especially with regards to the trailing `_PE1.fastq.gz` and `_PE2.fastq.gz`. 

We initialize six variables at the start:

`$directory`: location of the sequencing file

`$filename`: name of the gzipped FASTQ file

`$annotation_prefix`: Ensembl gene annotation GTF prefix; this is user-defined prefix above.

`$chrom_order_file`: a file with the FASTA reference entries of the genome ordered

`$UMI_length`: length of the UMI on the 5´ end of the paired end 1 read.

`$read_size`: read length minus UMI length.

`$cores`: number of cores for parallel processing

`$genome`: the genome FASTA file

`$genome_index`: the basename of the genome index files from `bowtie2-build`

`$prealign_rdna_index`: the basename of the prealign rDNA index files from `bowtie2-build`

```
directory=/Users/guertinlab/sequencing_run1 
filename=T47D_Starved_DMSO_rep1_PE1.fastq.gz
annotation_prefix=Homo_sapiens.GRCh38.104 
chrom_order_file=hg38.chrom.order.txt
UMI_length=8
read_size=62
cores=6
genome=hg38.fa
genome_index=hg38
prealign_rdna_index=human_rDNA
```

## Preprocessing 

Navigate to the folder that contains the sequencing file and gunzip the files.

```
cd $directory 
name=$(echo $filename | awk -F"_PE1.fastq.gz" '{print $1}')
echo $name
gunzip ${name}_PE*.fastq.gz
```

## Processing reads 

The first processing step is to remove adapter sequence and simultaneously discard reads that have insert sizes of one base. The 3´ adapter contains a UMI, which is sequenced prior to the adapter. Therefore, the vast majority of adapter/adapter ligation products have read lengths of exactly the length of the UMI. The option `-m $((UMI_length+2))` provides a one base buffer and discards reads with a length of the UMI + 1.

The fraction of reads that result from adapter/adapter ligation products is a useful metric to help determine the necessary read depth to achieve a certain aligned read depth. FASTQ files contain four lines per sequence entry, so we calculate this value by first counting the number of lines in the original FASTQ file using `wc -l` and divide that value by 4 using `awk '{print $1/4}'`. We perform the same operation on the output file of reads that were too short, in this case 0 or 1 base insertions. Finally we, divide the adapter/adapter ligation product value by the total and round to the hundredth with `$(echo "scale=2 ; $PE1_w_Adapter / $PE1_total" | bc)`.

This value varies widely depending upon whether a size selection was performed in the library preparation. We recently dropped the size selection step from the protocol in an effort to reduce bias against small RNA inserts [@sathyan2019improved]. If no size selection is performed, this value can be quite high; however, a high value does not indicate poor library quality. Because sequencing is relatively cheap, we tolerate up to 80% adapter/adapter ligation products. One needs to balance to cost of performing another experiment with sequencing uninformative adapter sequences. Later on we provide a formula for determining the required sequncing depth to result in a desired number of concordant aligned reads. So why the percentage of uninformative adapter reads is a useful number, if all other QC metrics point to high quality data, then we recommend further sequencing depth if less than 80% of the reads are adapter/adapter ligation products. 

Each quality control metric is distilled down to a single number that is printed to the file `${name}_QC_metrics.txt`. Alongside the value, we include the recommended threshold. We continue to print all metrics to this file and plot the data at the end of the workflow.

```
cutadapt --cores=$cores -m $((UMI_length+2)) -O 1 -a TGGAATTCTCGGGTGCCAAGG ${name}_PE1.fastq -o ${name}_PE1_noadap.fastq --too-short-output ${name}_PE1_short.fastq > ${name}_PE1_cutadapt.txt
cutadapt --cores=$cores -m $((UMI_length+10)) -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC ${name}_PE2.fastq -o ${name}_PE2_noadap.fastq --too-short-output ${name}_PE2_short.fastq > ${name}_PE2_cutadapt.txt

PE1_total=$(wc -l ${name}_PE1.fastq | awk '{print $1/4}')
PE1_w_Adapter=$(wc -l ${name}_PE1_short.fastq | awk '{print $1/4}')
AAligation=$(echo "scale=2 ; $PE1_w_Adapter / $PE1_total" | bc)

echo -e  "value\texperiment\tthreshold\tmetric" > ${name}_QC_metrics.txt
echo -e "$AAligation\t$name\t0.80\tAdapter/Adapter" >> ${name}_QC_metrics.txt
```

The next step removes reads with a length shorter than 10 bases and reverse complements the file so that the aligned read corresponds to the appropriate reference genome strand. This is necessary because the PRO-seq protocol sequences the 3´ end of the original nascent RNA as the paired end 1 read. 

```
seqtk seq -L $((UMI_length+10)) -r ${name}_PE1_noadap.fastq > ${name}_PE1_noadap_trimmed.fastq 
```

PRO-seq can have several independent reads that have the same genomic ends because promoter proximal pausing positions can be focused [@kwak2013precise] and the 5´end of the RNA is often the transcription start site. One cannot filter potential PCR duplicates based on whether two independent pairs of reads have identical paired end reads alignment. Therefore, we remove duplicate sequences from the FASTQ file based on the presence of the UMI. We effectively deduplicate the PE2 based on the presence of the PE1 UMI by pairing the PE1 and PE2 reads with `fastq_pair`. 

```
fqdedup -i ${name}_PE1_noadap_trimmed.fastq -o ${name}_PE1_dedup.fastq

#this variable is a near-optimal table size value for fastq_pair: 
PE1_noAdapter=$(wc -l ${name}_PE1_noadap.fastq | awk '{print $1/4}')

fastq_pair -t $PE1_noAdapter ${name}_PE1_noadap.fastq ${name}_PE2_noadap.fastq
```

## RNA integrity score

We measure RNA degradation by searching for overlap between paired end reads with `flash` and plotting the resultant histogram output with `insert_size.R`. We empirically found that there are fewer reads within the range of 10 - 20 than the range of 30 - 40 for high quality libraries [@smith2021peppro]. RNA only starts to protrude from the RNA Polymerase II exit channel at approximately 20 bases in length, so 20 bases of the nascent RNA is protected from degradation during the run on. Libraries with a substantial amount of degradation after the run on step are enriched for species in the range 10 - 20. A degradation ratio of less than 1 indicates a high quality library. Note that size selection to remove adapter/adapter ligation products will inflate this value because small RNAs are inevitably selected against when trying to remove only the adapter ligation band.    

```
flash -q --compress-prog=gzip --suffix=gz ${name}_PE1_noadap.fastq.paired.fq ${name}_PE2_noadap.fastq.paired.fq -o ${name}
insert_size.R ${name}.hist ${UMI_length}

rm ${name}_PE*_noadap.fastq.paired.fq
```

## Processing for alignment

The final processing step removes the UMI from both paired end reads and reverse complements the paired end 2 read.

```
seqtk trimfq -e ${UMI_length} ${name}_PE1_dedup.fastq > ${name}_PE1_processed.fastq
seqtk trimfq -e ${UMI_length} ${name}_PE2_noadap.fastq | seqtk seq -r - > ${name}_PE2_processed.fastq
```

## Remove reads aligning to rDNA
  
By first aligning to the rDNA, we can later estimate nascent RNA purity and avoid spurious read pile ups at region in the genome that non-uniquely align to both the rDNA locus and elsewhere in the genome. While between 70 - 80% of stable RNA is rRNA, generally less than 10% of the nascent RNA arises from rRNA. Even 10% of the library aligning to rDNA loci is extremely enriched, so any reads that map non-uniquely to both rDNA and non-rDNA regions in the genome result in huge artifactual spikes in the data if rDNA-aligned reads are not first removed. The `-f 0x4` flag of the `samtools fastq` command specifies that only unmapped reads should be included in the FASTQ output.    

```
bowtie2 -p $cores -x $prealign_rdna_index -U ${name}_PE1_processed.fastq 2>${name}_bowtie2_rDNA.log | samtools sort -n - | samtools fastq -f 0x4 - > ${name}_PE1.rDNA.fastq

#this removes PE2-aligned reads with a rDNA-aligned mate
reads=$(wc -l ${name}_PE1.rDNA.fastq | awk '{print $1/4}')
fastq_pair -t $reads ${name}_PE1.rDNA.fastq ${name}_PE2_processed.fastq
```

## Genome alignment
The last processing step for individual libraries is to align to the genome. Note that we reverse complemented both reads, so the `--rf` flag indicates that the reads are oriented opposite to typical RNA-seq data. The `samtools` commands convert the file to a compressed binary BAM format and sorts the reads. 

```
bowtie2 -p $cores --maxins 1000 -x $genome_index --rf -1 ${name}_PE1.rDNA.fastq.paired.fq -2 ${name}_PE2_processed.fastq.paired.fq 2>${name}_bowtie2.log | samtools view -b - | samtools sort - -o ${name}.bam
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
In order to calculate the rDNA alignment rate, we count the number of reads before and after removing rDNA aligned reads; their difference is the total number of rDNA-aligned reads. The command `samtools view -c -f 0x42` counts the concordant hg38-aligned paired end 1 reads. Lastly, we calculate the fraction of aligned reads that map to the rDNA genome and print it to the QC metrics file.    

```
PE1_prior_rDNA=$(wc -l ${name}_PE1_processed.fastq | awk '{print $1/4}')
PE1_post_rDNA=$(wc -l ${name}_PE1.rDNA.fastq | awk '{print $1/4}')
total_rDNA=$(echo "$(($PE1_prior_rDNA-$PE1_post_rDNA))") 

concordant_pe1=$(samtools view -c -f 0x42 ${name}.bam)
total_concordant=$(echo "$(($concordant_pe1+$total_rDNA))")

rDNA_alignment=$(echo "scale=2 ; $total_rDNA / $total_concordant" | bc)

echo -e "$rDNA_alignment\t$name\t0.10\trDNA Alignment Rate" >> ${name}_QC_metrics.txt
```

## Mappability rate
After all the processing, the vast majority of reads should map concordantly to the genome. Alignment rates for successful PRO-seq experiments are typically above 90%. Again, we count specific reads using `samtools` and by counts reads in FASTQ files, then calcuate alignment rate. We recommend the following site to help understand the meaning of samtools flags: https://broadinstitute.github.io/picard/explain-flags.html.

```
map_pe1=$(samtools view -c -f 0x40 -F 0x4 ${name}.bam)
pre_alignment=$(wc -l ${name}_PE1.rDNA.fastq.paired.fq | awk '{print $1/4}')
alignment_rate=$(echo "scale=2 ; $map_pe1 / $pre_alignment" | bc)

echo -e "$alignment_rate\t$name\t0.90\tAlignment Rate" >> ${name}_QC_metrics.txt
```

## Complexity and theoretical read depth

We developed `fqComplexity` to serve two purposes: 1) calculate the number of reads that are non-PCR duplicates as a metric for complexity; and 2) provide a formula and constants to calculate the theoretical read depth that will result in a user-defined number of concordant aligned reads. The equation accounts for all upstream processing steps. Over 10 milion concordantly aligned reads is typically sufficient if you have 3 or more replicates and a genome the size/gene density of the human genome. 

For the first metrics, the input FASTQ file is preprocessed with adapter/adapter products and inserts less than 10 removed; the UMI is still included in the FASTQ DNA sequence. The FASTQ file is subsampled into deciles and the intermediate files are deduplicated. The input and output numbers are logged in `${name}_complexity.log`. An asymptotic regression model is fit to the data and the total number of unique reads at 10 million read depth is printed on the resulting PDF plot. We recommend that at least 7.5 million reads are unique at a depth of 10 million (i.e. 75% of reads are unique if you were to align exactly 10 million reads).

```
fqComplexity -i ${name}_PE1_noadap_trimmed.fastq
```

We need two factors to derive the constants to calculate the theoretical read depth for a specified concordantly aligned read depth. First, we divide the total raw paired end 1 reads by processed paired end 1 reads that have adapter contamination and inserts less than 10 bases removed. The second value is the fraction of deduplicated reads that align concordantly to the non-rDNA genome.  

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

Invoke `fqComplexity` and specify the `-x` and `-y` options to fit an asymptotic regression model to the factor-scaled log file.

```
fqComplexity -i ${name}_PE1_noadap_trimmed.fastq -x $factorX -y $factorY
```

Counterintuitively, you can have a high quality and complex library that is not practical to sequence to further depth because the number of adapter/adapter
reads is too high. The QC metrics should be considered to determine whether the library is high quality. If the library is deemed high quality and low sequencing depth, use this equation to determine the practicality of increasing depth using the same libraries.   


## Get the reads in a BED
First we extract all the paired end 1 reads and separate reads based on their alignment to the plus or minus strand.

The software `seqOutBias` was originally developed to correct sequence bias from molecular genomics data. Although we are not correcting enzymatic sequence bias in this workflow, there are many features of `seqOutBias` that are useful. Note that we include the `--no-scale` option to not correct sequence bias. The software outputs a bigWig and BED6 file with strand information, but it also calculates mappability at the specified read length and excludes non-uniquely mappable reads. Lastly, invoking `--tail-edge` realigns the end of the read so that the exact position of RNA Polymerase is specified in the output BED and bigWig files.

```
samtools view -b -f 0x40 ${name}.bam > ${name}_PE1.bam
samtools view -bh -F 0x14 ${name}_PE1.bam > ${name}_PE1_plus.bam
samtools view -bh -f 0x10 ${name}_PE1.bam > ${name}_PE1_minus.bam
    
seqOutBias $genome ${name}_PE1_plus.bam --no-scale --bed ${name}_PE1_plus.bed --bw=${name}_PE1_plus.bigWig --tail-edge --read-size=$read_size --stranded --bed-stranded-positive
seqOutBias $genome ${name}_PE1_minus.bam --no-scale --bed ${name}_PE1_minus.bed --bw=${name}_PE1_minus.bigWig --tail-edge --read-size=$read_size --stranded --bed-stranded-positive

cat ${name}_PE1_plus.bed ${name}_PE1_minus.bed | sort -k1,1 -k2,2n > ${name}_PE1_signal.bed
```


## Run on efficiency

RNA polymerases that are associated with gene bodies efficiently incorporate nucleotides during the run on reaction under most conditions, but promoter proximal paused RNA polymerase require high salt or detergent to run on efficiently [@rougvie1988rna; @core2012defining]. Therefore, the pause index is used to quantify run on efficiency. Pause index is the density of signal in the promoter-proximal pause region divided by density in the gene body. However, since pause windows are user-defined and variable, pause indices can differ substantially between metrics. There are many exon 1 gene annotations depending on gene isoforms and the upstream most annotated TSS is not necessarily the prominently transcribed isoform. It is common practice to choose the upstream most TSS, but this will cause the pause index to be artificially deflated. Here, we define the pause window for a gene as position 20 - 120 downstream of the most prominent TSS. The most prominent TSS is determined by calculating the density in this 20 - 120 window for all annotated TSSs for each gene and choosing the TSS upstream of the most RNA-polymerase dense region for each gene.   

```
coverageBed -sorted -counts -s -a $annotation_prefix.pause.bed -b ${name}_PE1_signal.bed -g $chrom_order_file | awk '$7>0' | sort -k5,5 -k7,7nr | sort -k5,5 -u > ${name}_pause.bed

#discard anything with chr and strand inconsistencies
join -1 5 -2 5 ${name}_pause.bed $annotation_prefix.bed | awk '{OFS="\t";} $2==$8 && $6==$12 {print $2, $3, $4, $1, $6, $7, $9, $10}' | awk '{OFS="\t";} $5 == "+" {print $1,$2+480,$8,$4,$6,$5} $5 == "-" {print $1,$7,$2 - 380,$4,$6,$5}' |  awk  '{OFS="\t";} $3>$2 {print $1,$2,$3,$4,$5,$6}' | sort -k1,1 -k2,2n > ${name}_pause_counts_body_coordinates.bed

#column ten is Pause index
coverageBed -sorted -counts -s -a ${name}_pause_counts_body_coordinates.bed -b ${name}_PE1_signal.bed -g $chrom_order_file | awk '$7>0' | awk '{OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,$5/100,$7/($3 - $2)}' | awk '{OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$8/$9}' > ${name}_pause_body.bed

pause_index.R ${name}_pause_body.bed
```

## Estimate nascent RNA purity with exon / intron density ratio

RNA-seq primarily measures mature transcripts, so exons density far exceeds intron density. However, these densities are comparable within gene bodies for nascent RNA-seq. For calculation of this metric, the first exon is excluded because pausing occurs in this region and artifically inflates the exon density. Calculating the exon density to intron density ratio complements rDNA alignment rate as a metric to quantify nascent RNA purity. 

```
coverageBed -sorted -counts -s -a $annotation_prefix.introns.bed -b ${name}_PE1_signal.bed -g $chrom_order_file  | awk '$7>0' | awk '{OFS="\t";} {print $1,$2,$3,$5,$5,$6,$7,($3 - $2)}' > ${name}_intron_counts.bed

coverageBed -sorted -counts -s -a $annotation_prefix.no.first.exons.named.bed -b ${name}_PE1_signal.bed -g $chrom_order_file | awk '$7>0' | awk '{OFS="\t";} {print $1,$2,$3,$4,$4,$6,$7,($3 - $2)}' > ${name}_exon_counts.bed

exon_intron_ratio.R ${name}_exon_counts.bed ${name}_intron_counts.bed
```
## Remove intermediate files and zip raw sequencing files

Calculating these quality control metrics necessitates many intermediate files. Many files are unused output from various processing steps, but many are only used briefly for counting reads before and after processing. FASTQ files are large and rarely used in downstream analyses, so the following code chunk removes intermediate FASTQ files and compresses the original files. 
```
rm ${name}_PE1_short.fastq
rm ${name}_PE2_short.fastq
rm ${name}_PE1_noadap.fastq
rm ${name}_PE2_noadap.fastq
rm ${name}_PE1_noadap_trimmed.fastq
rm ${name}_PE1_dedup.fastq
rm ${name}_PE1_processed.fastq
rm ${name}_PE2_processed.fastq
rm ${name}_PE1.rDNA.fastq.paired.fq
rm ${name}_PE1.rDNA.fastq.single.fq
rm ${name}_PE2_processed.fastq.paired.fq
rm ${name}_PE2_processed.fastq.single.fq
rm ${name}_PE1_noadap.fastq.single.fq
rm ${name}_PE2_noadap.fastq.single.fq
rm ${name}.extendedFrags.fastq.gz
rm ${name}.notCombined_2.fastq.gz
rm ${name}.notCombined_1.fastq.gz
gzip ${name}_PE1.fastq
gzip ${name}_PE2.fastq
```

## Process all files in series

It is a useful exercise to run through the code chunks above individually and look at each output to gain further understanding of each step. A full understanding of the workflow allows the end user to modify the steps to account for modifications in the PRO-seq protocol. However, automation of routine processing and analysis is more practical once a workflow is established. Below, we provide a simple shell script loop that will process each set of files in series. This loop can be easily adapted to perform all processing in parallel using a job scheduler and submission of a batch script for each set of input files.



```
#initialize variables
directory=/Users/guertinlab/sequencing_run1_series 
annotation_prefix=Homo_sapiens.GRCh38.104 
chrom_order_file=hg38.chrom.order.txt
UMI_length=8
read_size=62
cores=6
genome=hg38.fa
genome_index=hg38
prealign_rdna=human_rDNA

cd $directory 

for filename in *PE1.fastq.gz
do
    name=$(echo $filename | awk -F"_PE1.fastq.gz" '{print $1}')
    echo $name
    echo 'unzipping raw' $name 'files'
    gunzip ${name}_PE*.fastq.gz
    echo 'removing dual adapter ligations and calculating the fraction of adapter/adapters in' $name
    cutadapt --cores=$cores -m $((UMI_length+2)) -O 1 -a TGGAATTCTCGGGTGCCAAGG ${name}_PE1.fastq -o ${name}_PE1_noadap.fastq --too-short-output ${name}_PE1_short.fastq > ${name}_PE1_cutadapt.txt
    cutadapt --cores=$cores -m $((UMI_length+10)) -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC ${name}_PE2.fastq -o ${name}_PE2_noadap.fastq --too-short-output ${name}_PE2_short.fastq > ${name}_PE2_cutadapt.txt
    PE1_total=$(wc -l ${name}_PE1.fastq | awk '{print $1/4}')
    PE1_w_Adapter=$(wc -l ${name}_PE1_short.fastq | awk '{print $1/4}')
    AAligation=$(echo "scale=2 ; $PE1_w_Adapter / $PE1_total" | bc)
    echo -e  "value\texperiment\tthreshold\tmetric" > ${name}_QC_metrics.txt
    echo -e "$AAligation\t$name\t0.80\tAdapter/Adapter" >> ${name}_QC_metrics.txt
    echo 'removing short RNA insertions and reverse complementing in' $name
    seqtk seq -L $((UMI_length+10)) -r ${name}_PE1_noadap.fastq > ${name}_PE1_noadap_trimmed.fastq
    echo 'removing PCR duplicates from' $name
    fqdedup -i ${name}_PE1_noadap_trimmed.fastq -o ${name}_PE1_dedup.fastq
    PE1_noAdapter=$(wc -l ${name}_PE1_noadap.fastq | awk '{print $1/4}')
    fastq_pair -t $PE1_noAdapter ${name}_PE1_noadap.fastq ${name}_PE2_noadap.fastq
    echo 'calculating and plotting RNA insert sizes from' $name
    flash -q --compress-prog=gzip --suffix=gz ${name}_PE1_noadap.fastq.paired.fq ${name}_PE2_noadap.fastq.paired.fq -o ${name}
    insert_size.R ${name}.hist ${UMI_length}
    echo 'trimming off the UMI from' $name
    seqtk trimfq -e ${UMI_length} ${name}_PE1_dedup.fastq > ${name}_PE1_processed.fastq
    seqtk trimfq -e ${UMI_length} ${name}_PE2_noadap.fastq | seqtk seq -r - > ${name}_PE2_processed.fastq
    echo 'aligning' $name 'to rDNA and removing aligned reads'
    bowtie2 -p $cores -x $prealign_rdna -U ${name}_PE1_processed.fastq 2>${name}_bowtie2_rDNA.log | samtools sort -n - | samtools fastq -f 0x4 - > ${name}_PE1.rDNA.fastq
    reads=$(wc -l ${name}_PE1.rDNA.fastq | awk '{print $1/4}')
    fastq_pair -t $reads ${name}_PE1.rDNA.fastq ${name}_PE2_processed.fastq
    echo 'aligning' $name 'to the genome'
    bowtie2 -p $cores --maxins 1000 -x $genome_index --rf -1 ${name}_PE1.rDNA.fastq.paired.fq -2 ${name}_PE2_processed.fastq.paired.fq 2>${name}_bowtie2.log | samtools view -b - | samtools sort - -o ${name}.bam
    PE1_prior_rDNA=$(wc -l ${name}_PE1_processed.fastq | awk '{print $1/4}')
    PE1_post_rDNA=$(wc -l ${name}_PE1.rDNA.fastq | awk '{print $1/4}')
    total_rDNA=$(echo "$(($PE1_prior_rDNA-$PE1_post_rDNA))") 
    echo 'calculating rDNA and genomic alignment rates for' $name
    concordant_pe1=$(samtools view -c -f 0x42 ${name}.bam)
    total_concordant=$(echo "$(($concordant_pe1+$total_rDNA))")
    rDNA_alignment=$(echo "scale=2 ; $total_rDNA / $total_concordant" | bc)
    echo -e "$rDNA_alignment\t$name\t0.20\trDNA Alignment Rate" >> ${name}_QC_metrics.txt
    map_pe1=$(samtools view -c -f 0x40 -F 0x4 ${name}.bam)
    pre_alignment=$(wc -l ${name}_PE1.rDNA.fastq.paired.fq | awk '{print $1/4}')
    alignment_rate=$(echo "scale=2 ; $map_pe1 / $pre_alignment" | bc)
    echo -e "$alignment_rate\t$name\t0.90\tAlignment Rate" >> ${name}_QC_metrics.txt
    echo 'plotting and calculating complexity for' $name
    fqComplexity -i ${name}_PE1_noadap_trimmed.fastq
    echo 'calculating and plotting theoretical sequencing depth to achieve a defined number of concordantly aligned reads for' $name
    PE1_total=$(wc -l ${name}_PE1.fastq | awk '{print $1/4}')
    PE1_noadap_trimmed=$(wc -l ${name}_PE1_noadap_trimmed.fastq | awk '{print $1/4}')
    factorX=$(echo "scale=2 ; $PE1_total / $PE1_noadap_trimmed" | bc)
    echo fraction of reads that are not adapter/adapter ligation products or below 10 base inserts
    echo $factorX | awk '{print 1/$1}'
    PE1_dedup=$(wc -l ${name}_PE1_dedup.fastq | awk '{print $1/4}')
    factorY=$(echo "scale=2 ; $concordant_pe1 / $PE1_dedup" | bc)
    fqComplexity -i ${name}_PE1_noadap_trimmed.fastq -x $factorX -y $factorY
    echo 'Separating paired end reads and creating genomic BED and bigWig intensity files for' $name
    samtools view -b -f 0x40 ${name}.bam > ${name}_PE1.bam
    samtools view -bh -F 0x14 ${name}_PE1.bam > ${name}_PE1_plus.bam
    samtools view -bh -f 0x10 ${name}_PE1.bam > ${name}_PE1_minus.bam
    seqOutBias $genome ${name}_PE1_plus.bam --no-scale --bed ${name}_PE1_plus.bed --bw=${name}_PE1_plus.bigWig --tail-edge --read-size=$read_size --stranded --bed-stranded-positive
    seqOutBias $genome ${name}_PE1_minus.bam --no-scale --bed ${name}_PE1_minus.bed --bw=${name}_PE1_minus.bigWig --tail-edge --read-size=$read_size --stranded --bed-stranded-positive
    cat ${name}_PE1_plus.bed ${name}_PE1_minus.bed | sort -k1,1 -k2,2n > ${name}_PE1_signal.bed
    echo 'Calculating pause indices for' $name  
    coverageBed -sorted -counts -s -a $annotation_prefix.pause.bed -b ${name}_PE1_signal.bed -g $chrom_order_file | awk '$7>0' | sort -k5,5 -k7,7nr | sort -k5,5 -u > ${name}_pause.bed
    join -1 5 -2 5 ${name}_pause.bed $annotation_prefix.bed | awk '{OFS="\t";} $2==$8 && $6==$12 {print $2, $3, $4, $1, $6, $7, $9, $10}' | awk '{OFS="\t";} $5 == "+" {print $1,$2+480,$8,$4,$6,$5} $5 == "-" {print $1,$7,$2 - 380,$4,$6,$5}' |  awk  '{OFS="\t";} $3>$2 {print $1,$2,$3,$4,$5,$6}' | sort -k1,1 -k2,2n > ${name}_pause_counts_body_coordinates.bed
    coverageBed -sorted -counts -s -a ${name}_pause_counts_body_coordinates.bed -b ${name}_PE1_signal.bed -g $chrom_order_file | awk '$7>0' | awk '{OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,$5/100,$7/($3 - $2)}' | awk '{OFS="\t";} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$8/$9}' > ${name}_pause_body.bed
    pause_index.R ${name}_pause_body.bed
    echo 'Calculating exon density / intron density as a metric for nascent RNA purity for' $name
    coverageBed -sorted -counts -s -a $annotation_prefix.introns.bed -b ${name}_PE1_signal.bed -g $chrom_order_file | awk '$7>0' | awk '{OFS="\t";} {print $1,$2,$3,$5,$5,$6,$7,($3 - $2)}' > ${name}_intron_counts.bed
    coverageBed -sorted -counts -s -a $annotation_prefix.no.first.exons.named.bed -b ${name}_PE1_signal.bed -g $chrom_order_file | awk '$7>0' | awk '{OFS="\t";} {print $1,$2,$3,$4,$4,$6,$7,($3 - $2)}' > ${name}_exon_counts.bed
    exon_intron_ratio.R ${name}_exon_counts.bed ${name}_intron_counts.bed
    #clean up intermediate files and gzip
    rm ${name}_PE1_short.fastq
    rm ${name}_PE2_short.fastq
    rm ${name}_PE1_noadap.fastq
    rm ${name}_PE2_noadap.fastq
    rm ${name}_PE1_noadap_trimmed.fastq
    rm ${name}_PE1_dedup.fastq
    rm ${name}_PE1_processed.fastq
    rm ${name}_PE2_processed.fastq
    rm ${name}_PE1_noadap.fastq.paired.fq   
    rm ${name}_PE2_noadap.fastq.paired.fq
    rm ${name}_PE1.rDNA.fastq.paired.fq
    rm ${name}_PE1.rDNA.fastq.single.fq
    rm ${name}_PE2_processed.fastq.paired.fq
    rm ${name}_PE2_processed.fastq.single.fq
    rm ${name}_PE1_noadap.fastq.single.fq
    rm ${name}_PE2_noadap.fastq.single.fq
    rm ${name}.extendedFrags.fastq.gz
    rm ${name}.notCombined_2.fastq.gz
    rm ${name}.notCombined_1.fastq.gz
    gzip ${name}_PE1.fastq
    gzip ${name}_PE2.fastq
done

```
## Plot all QC metrics

Individual detailed quality controls plots provide valuable information about the data, but each metric can be distilled to a single value. We empirically determined thresholds for each value that constitute acceptable libraries. These thresholds are not considered absolute and should only be used as guidelines. Below, we concatenate all the QC metrics for the experiments and plot the results. Thresholds lines are included on the plots. One can quickly glance at the plot to determine whether the QC values fall within the light green shaded area. If values are outside the light green shaded area, then we recommend looking back at the more detailed QC plots to diagnose possible issues with the libraries.  

```
cat *_QC_metrics.txt | awk '!x[$0]++' > project_QC_metrics.txt 

plot_all_metrics.R project_QC_metrics.txt Estrogen_treatment_PRO
```

## Differential expression with DESeq2

Differential expression analysis is a common first step after routine RNA-seq and PRO-seq data processing. Below we outline the `bedtools` commands to count reads within gene annotations and we provide an `R` script for differentially expression analysis with `DESeq2`. The script also plots the fold change between conditions and mean expression level for each gene.  For simplicity we use the most upstream transcription start site and most downstream transcription termination site for annotations, but we recommend more sophisticated methods to define primary transcripts [@anderson2020defining; zhao2021deconvolution]. The `R` script requires three ordered arguments: 1) a file with the signal counts for each gene in every even row, 2) the prefix for the baseline experimental condition for which to compare (often termed "untreated"), 3) prefix name for the output PDF plot. 


```
annotation_prefix=Homo_sapiens.GRCh38.104 
chrom_order_file=hg38.chrom.order.txt
sort -k1,1 -k2,2n $annotation_prefix.bed > $annotation_prefix.sorted.bed 

for filename in *_PE1_signal.bed
do
    name=$(echo $filename | awk -F"_PE1_signal.bed" '{print $1}')
    echo -e  "\t${name}" > ${name}_gene_counts.txt
    coverageBed -sorted -counts -s -a $annotation_prefix.sorted.bed  -b $filename -g $chrom_order_file | awk '{OFS="\t";} {print $4,$7}' >> ${name}_gene_counts.txt
done
paste -d'\t' *_gene_counts.txt > Estrogen_treatment_PRO_gene_counts.txt

differential_expression.R Estrogen_treatment_PRO_gene_counts.txt T47D_Starved_DMSO Estrogen_treatment

```
## Conclusions

Analyses are presented in a deconstructed manner to provide flexibility to researchers who wish to develop their own workflows and pipelines, or as Captain Barbossa succinctly stated: “The code is more of what you call _guidelines_ than actual rules.”
[@martins2018universal]


# References
