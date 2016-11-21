////////////////////////////////////////////////////////////
// GATK-based 'best' practices variant-calling pipeline for 
// Plasmodium falciparum.
// 
// This pipeline is a port of the VLSCI whole genome variant
// calling pipeline  for from this Git repo to Bpipe:
// 
//     https://github.com/claresloggett/variant_calling_pipeline/
// And modified from
//     https://github.com/ssadedin/variant_calling_pipeline/
//
//
// This pipeline attempts to recreate the WGS processing pipeline developed
// by the Malaria Genomics Consortium for Plasmodium falciparum. 
//
////////////////////////////////////////////////////////////

// Create this file by copying config.groovy.template and editing
load 'config.groovy'

// All the core pipeline stages in the pipeline
load 'pipeline_stages_config.groovy'


run {  
    "%_*_R*.fastq.gz" * [   makeUBAM  +
                  fastqc_unmapped +   
                  markAdapters +
                  remapBWA + 
                  realignIntervals + 
                  realignIndels + 
                  dedup + 
                  bqsrPass1 + 
                  bqsrPass2 + 
                  bqsrApply + 
                  bqsrCheck + 
                  fastqc_mapped + 
                  alignment_metrics +
                  insert_metrics + 
                  genomecoverage +
                  perbasecoveragepf +
                  callVariants ] + 
    // after each BAM processed jointly genotype then call biallelic variants
    [
                  combineGVCF + 
                  vqsrGenerate + 
                  vqsrApply + 
                  annotate +
                  regions +
                  barcode + 
                  keepSNPs + 
                  filterSNPs +
                  cleanVCF +
                  extractAnno + 
                  indexVCF

     ] 
}
