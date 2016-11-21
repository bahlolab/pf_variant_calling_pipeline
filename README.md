# Plasmodium falciparum WGS variant calling pipeline
This is my attempt to replicate the pipeline used by the 
[Malaria Genomics Consortium](https://www.malariagen.net/projects/parasite/pf) 
for calling single nucleotide variants (SNVs) from paired-end whole-genome 
sequencing data.

It is a [bpipe](http://docs.bpipe.org/Overview/Introduction/) pipeline adapted 
from the [example](https://github.com/ssadedin/variant_calling_pipeline/)
[VLSCI](https://github.com/claresloggett/variant_calling_pipeline/) pipeline for
calling variants from WGS data in humans. 

## Overview

The Malaria Gen pipeline follows the GATK best practices guideline and takes
as paired FASTQ files. 

It extracts the fastq files by read-group and marks adapter sequences. Following
the reads for each sample are mapped using bwa mem -M and remerged into one BAM file.
Then the usual GATK steps: indel relalignment, duplicate marking, base quality score
recalibration are performed. QC steps such as coverage analysis are performed.
Variant calling is down for each sample using
GATK Haplotype Caller in discovery mode and outputting gVCF files.

The gVCF files are collectively genotyped and variant-quality score recalibration
is applied. The VCF is then annotated using snpEff and custom annotations are applied
using bcftools. The final VCF consists of biallelic SNPs that are contained in 
core regions and have VQSLOD > 0. You can alter these settings if you require less
sensitivity. 

### What's included?
We have included the _P.falciparum_ reference and genome annotations files from
MalariaGen. From these files we have constructed the bwa indexes + snpEff databases
for mapping and annotation, respecitvely. 

For variant calling and annotation, we have included the vcf files from the _P.falciparum_
genetic crosses community project to perform base-quality score recalibration and 
variant-quality score recalibration. We have also used data from this project
to annotate SNPs in hypervariable or subtelomeric regions. We also constructed
a list of known drug resistance genes to assess and positions for geographically
informative SNPs from [Neafsey et. al, 2008](http://www.ncbi.nlm.nih.gov/pubmed/19077304).

### What's outputted?
This pipeline constructs an 'analysis' ready fully annotated VCF file based
on cleaning recommendations implemented by the Malaria Genomics Consortium. 
This VCF file contains custom annotations based on wheter a SNP is party of the
Neafsey barcode. We also extract a genomic data storage file for manipulating
VCF files in R. 

The pipeline is modular so if extra steps are required for your project
it is quite easy to implement them as additional grrovy variable.


## Requirements and usage
There are a number of software requirements, which you should ensure are 
satisfied before running the pipeline. These need to be configured 
in the file called 'config.groovy'.

At the moment we use the following software in our pipeline

* bwa mem 0.7.10-r789
* PicardTools 1.99
* FastQC  
* GATK v3.5
* SnpEff 4.11 
* bcftools 1.3
* R 3.2.4 (with SeqArray library installed)

To use the pipeline first install bpipe (version 0.9.8.5 or later)
and then simply enter 

```{bash}
bpipe run -n 4 pipeline.groovy sample_01_R1.fq.gz sample_02_R2.fq.gz
```
to run the the pipeline concurrently using four cores. We have
also found that for some computationally intensive steps it's worth
limiting the availble memory to all processes by reducing MAX_JAVA_MEM.

```{bash}
MAX_JAVA_MEM=4g bpipe run -n 4 pipeline.groovy *.bam
```

# Appendix

All the genetic crosses, fasta and gff files for this pipeline were obtained from
Malaria Gen Pf3k v5 pilot data.
PFalciparum - ftp://ngs.sanger.ac.uk/production/pf3k/release_5/
