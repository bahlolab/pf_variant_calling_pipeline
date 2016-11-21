////////////////////////////////////////////////////////
//
// Pipeline stage definitions for example WGS variant calling 
// pipeline. See pipeline.groovy for more information.
// 
////////////////////////////////////////////////////////
createBWAIndex = {
    doc "Run BWA index on fasta if not available."

    exec """
        bwa index -a bwtsw -p Pfalciparum.genome $REF
        mv Pfalciparum.genome.* $REFBASE/fasta/

    """

}

fastaIndex = {
    doc "index fasta file if not available"
    exec "samtools faidx $REF"
}

createDictionary = {
    doc "Create Sequence dictionary"
    output.dir="$REFBASE/fasta"
    from("fasta") to("dict") {
        exec "java -jar $PICARD_HOME/CreateSequenceDictionary.jar R=$REF O=$output.dict"
    }
}

makeUBAM = {
    doc "Create uBAM from input fastq file using revert sam. Note we clean bam (revert mapQ=0 for unmapped reads) file before reverting."
    output.dir="$REFBASE/bam"
    branch.sample="$branch.name"
    produce("$sample" + ".raw.bam") {
        exec """            
            $PICARD_HOME/FastqToSam \
                F1=$input1.gz \
                F2=$input2.gz \
                O=$output \
                SM=$sample \
                RG=$sample \
                PU=\$(zcat $input1.gz | head -n1 | cut -d":" -f3,6,7 | sed 's/:/./g' | sed 's/\\s1//g') \
                PLATFORM=ILLUMINA \
                SEQUENCING_CENTER=WEHI \
                SORT_ORDER=queryname \
                TMP_DIR=$TMPDIR 
        """
    }
}

fastqc_unmapped = {
    doc "Run FASTQC to generate QC metrics for the reads post-alignment"
    output.dir = "$REFBASE/qc"
    transform('.bam')  to('_fastqc.zip')  {
        exec "fastqc --quiet -o ${output.dir} -f bam $input.bam"
    }
}

@filter("markilluminaadapters")
markAdapters = {
    doc "Mark Illumina Adapters with PicardTools"
    output.dir="$REFBASE/ubam"
    def adapter_file="$output.dir" + "/" + "$sample" +"adapters.txt"
    uses(GB:8) {
        exec """ 
            $PICARD_HOME/MarkIlluminaAdapters \
                I=$input.bam \
                O=$output.bam \
                PE=true \
                M=$adapter_file \
                TMP_DIR=$TMPDIR

        """
    }
}

@preserve
remapBWA = {
    doc "Create merged BAM alignment from unmapped BAM file"
    output.dir="$REFBASE/bam"
    def ubam="$REFBASE/bam/"+"$sample"+".raw.bam"
    uses(threads:4,GB:16) {
    exec """
        $PICARD_HOME/SamToFastq I=$input.bam FASTQ=/dev/stdout \
            CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
            TMP_DIR=$TMPDIR | \
        bwa mem -M -t 4 -p $INDEXES /dev/stdin | \
        $PICARD_HOME/MergeBamAlignment \
            ALIGNED_BAM=/dev/stdin \
            UNMAPPED_BAM=$ubam\
            OUTPUT=$output.bam \
            R=$REF CREATE_INDEX=true ADD_MATE_CIGAR=true \
            CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
            INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
            PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
            TMP_DIR=$TMPDIR


        """
    }

}

@transform(".intervals")
realignIntervals = {
    output.dir="$REFBASE/bam"
    exec """
        java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
          -T RealignerTargetCreator
          -R $REF 
          -I $input.bam
          -log $LOG
          -o $output.intervals
    """
}

@filter("realignindels")
realignIndels = {
    output.dir="$REFBASE/bam"
    exec """
        java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar 
          -T IndelRealigner 
          -R $REF 
          -I $input.bam
          -targetIntervals $input.intervals 
          -log $LOG
          -o $output.bam
    """
}

@filter("dedup")
dedup = {
    doc "Mark Duplicates with PicardTools"
    output.dir="$REFBASE/bam"
    def metrics="$REFBASE/qc/" + "$sample" +".dup.metrics.txt" 
    exec """
        $PICARD_HOME/MarkDuplicates
             INPUT=$input.bam
             REMOVE_DUPLICATES=true 
             VALIDATION_STRINGENCY=LENIENT 
             AS=true 
             METRICS_FILE=$metrics
             CREATE_INDEX=true 
             OUTPUT=$output.bam
    """
}

bqsrPass1 = {
    doc "Recalibrate base qualities in a BAM file so that quality metrics match actual observed error rates"
    output.dir="$REFBASE/bam"
    exec """
            java -Xmx6g -jar $GATK/GenomeAnalysisTK.jar
                -T BaseRecalibrator 
                -R $REF  
                -I $input.bam  
                -knownSites $REFSNP1  
                -knownSites $REFSNP2 
                -knownSites $REFSNP3
                -o $output.pass1.table 
        """
    
}


bqsrPass2 = {
    doc "Recalibrate base qualities in a BAM file so that quality metrics match actual observed error rates"
    output.dir="$REFBASE/bam"
    exec """
            java -Xmx6g -jar $GATK/GenomeAnalysisTK.jar 
                -T BaseRecalibrator 
                -R $REF 
                -I $input.bam 
                -knownSites $REFSNP1  
                -knownSites $REFSNP2 
                -knownSites $REFSNP3 
                -BQSR $input.pass1.table
                -o $output.pass2.table 
        """
}

@preserve
bqsrApply = {
    doc "Apply BQSR to input BAM file."
    output.dir="$REFBASE/bam"
    exec """
        java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar  
            -T PrintReads  
            -R $REF 
            -I $input.bam   
            -BQSR $input.pass1.table  
            -o $output.bam
    """ 
}

bqsrCheck = {
    doc "Compare pre and post base-quality scores from recalibration"
    output.dir="$REFBASE/qc"
    transform('.pass1.table', '.pass2.table') to('.csv') {
        exec """
            java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
              -T AnalyzeCovariates
              -R $REF
              -before $input.pass1.table 
              -after $input.pass2.table 
              -csv $output
        """
    }
}


// QC metrics at BAM level
fastqc_mapped = {
    doc "Run FASTQC to generate QC metrics for the reads post-alignment"
    output.dir = "$REFBASE/qc"
    transform('.bam')  to('_fastqc.zip')  {
        exec "fastqc --quiet -o ${output.dir} -f bam_mapped $input.bam"
    }
}

@transform(".alignment_metrics")
alignment_metrics = {
    doc "Collect alignment summary statistics"
    output.dir="qc"
    produce("$sample" + ".alignment_metrics") {
        exec """ 
            $PICARD_HOME/CollectAlignmentSummaryMetrics \
                R=$REF \
                I=$input.bam \
                O=$output
        """
    }
}

@transform(".insert_metrics")
insert_metrics = {
    doc "Collect insert size metrics"
    output.dir="qc"
    def histogram="$output.dir" + "/" + "$sample" + "is_distribution.pdf"
    produce("$sample" + ".insert_metrics") {
        exec """
            $PICARD_HOME/CollectInsertSizeMetrics \
                I=$input.bam \
                O=$output.insert_metrics \
                H=$histogram 
        """
    }

}

// Depth of Coverage Metrics
genomecoverage = {
    output.dir = "qc"
    produce("$sample" + ".bw") {
        exec "bamCoverage -b $input.bam --minMappingQuality 30 --samFlagInclude 64 -o $output"
    }
} 

perbasecoveragepf = {
    output.dir = "qc"
    produce("$sample" + ".pf_coverage") {
            exec """
                bedtools coverage -abam $input.bam -b data/fasta/pf.contigs.sorted.bed -hist > $output
            """
    }
}

// Variant calling
callVariants = {
    doc "Call SNPs/SNVs using GATK Haplotype Caller, produces .g.vcf"
    output.dir="$REFBASE/variants"
    transform(".bam") to(".g.vcf") {
        exec """
            java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar
                -T HaplotypeCaller 
                -R $REF
                -I $input.bam
                --emitRefConfidence GVCF
                -gt_mode DISCOVERY
                -o $output

        """
    }
}

@preserve
combineGVCF = {
    doc "Jointly genotype gVCF files"
    output.dir="$REFBASE/variants_combined"
    def gvcfs = "--variant " + "$inputs.g.vcf".split(" ").join(" --variant ")
    produce("combined_variants.vcf") {
        exec """
            java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar
                -T GenotypeGVCFs
                -R $REF
                $gvcfs
                -o $output
    """
    }
}

// variant recalibration
vqsrGenerate = {
    doc "Generate variant quality score recalibration. Note requires GATK version 3.5"
    output.dir="$REFBASE/variants_combined"
    exec """
         java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar 
            -T VariantRecalibrator 
            -R $REF 
            -input $input.vcf
            -resource:cross1,training=true,truth=true,known=false $REFSNP1
            -resource:cross2,training=true,truth=true,known=false $REFSNP2
            -resource:cross3,training=true,truth=true,known=false $REFSNP3
            -an QD
            -an MQ
            -an FS 
            -an SOR 
            -an DP
            -mode SNP 
            -mG 8
            -MQCap 70 
            -recalFile $output.recal 
            -tranchesFile $output.tranches 
            -rscriptFile $output.plots.R

    """
}

vqsrApply = {
    doc "Apply variant quality score recalibration"
    output.dir="$REFBASE/variants_combined"
    exec """
        java -jar $GATK/GenomeAnalysisTK.jar 
            -T ApplyRecalibration 
            -R $REF  
            -input $input.vcf  
            -mode SNP  
            --ts_filter_level 99.0  
            -recalFile $input.recal  
            -tranchesFile $input.tranches  
            -o $output.vcf
    """
}

// variant annotation using snpEff
annotate = {
    doc "Annotate variants using snpEff"
    output.dir="$REFBASE/variants_combined"
    exec """ 
        java -Xmx8g -jar  $SNPEFF_HOME/snpEff.jar
            -c $SNPEFF_CONFIG
            -no-downstream 
            -no-upstream 
            Pf3D7v3
            $input.vcf > $output.vcf ;  

    """
}

regions = {
    doc "Include core regions from Pf genetic crosses version 1"
    output.dir="$REFBASE/variants_combined"

    exec """

        bcftools annotate -a $CORE_REGIONS -h $CORE_REGIONS_HDR -Ov -o $output.vcf -c CHROM,FROM,TO,RegionType $input.vcf


    """
}

barcode = {
    doc "Annotate global barcode SNPs from Neafsey et al., 2008"
    output.dir="$REFBASE/variants_combined"

    exec """
        bcftools annotate -a $BARCODE -h $BARCODE_HDR -Ov -o $output.vcf -c CHROM,FROM,TO,GlobalBarcode $input.vcf

    """ 
}

// select only biallelic SNPs
keepSNPs= {
    doc "Kepp only SNPs in VCF file"
    output.dir="$REFBASE/variants_combined"
    produce("biallelic_only.vcf") {
        exec """
            java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
                -T SelectVariants \
                -R $REF
                -V $input.vcf \
                -o $output \
                -selectType SNP
                --restrictAllelesTo BIALLELIC
    """
    }
}

// variant filtering annotations
@filter("filter")
filterSNPs = {
    doc "Annotate VCF file with additional filters at the variant level"
    output.dir="$REFBASE/variants_combined"
    exec """
        java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
            -T VariantFiltration 
            -R $REF
            --filterName LowQualVQ -filter "VQSLOD <= 0.0"
            --filterName NotCore -filter  "RegionType != 'Core'"
            --variant $input.vcf
            -log $LOG 
            -o $output.vcf
    """
}

// apply GATK SelectVariants to filter Low Qual regions
@preserve 
cleanVCF = {
    doc "Clean VCF for analysis ready."
    output.dir="cache"
    produce("final_snps.vcf") {
        exec """
            java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar
                -T SelectVariants
                -R $REF
                --variant $input.vcf
                -o $output
                -select 'vc.isNotFiltered()'
        """
    }
}


extractAnno = {
    doc "Use snpSift to extract annotations as plain text file"
    output.dir="cache"
    produce("final_annotations.txt") {
        exec """
            cat $input.vcf | perl $SNPEFF_HOME/scripts/vcfEffOnePerLine.pl |\
            java -Xmx4g -jar $SNPEFF_HOME/SnpSift.jar extractFields -e "NA" - \
                CHROM POS REF ALT "ANN[*].ALLELE" GlobalBarcode DP "ANN[*].GENEID" "ANN[*].BIOTYPE" "ANN[*].EFFECT"  "ANN[*].HGVS_P" "ANN[*].ERRORS" > $output
        """

    }
}

indexVCF = {
    exec "./vcftools_prepare.sh $input.vcf"
}


