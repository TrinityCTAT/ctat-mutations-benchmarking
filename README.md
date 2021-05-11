### Coversion of VCF output file from variant calling pipeline to standard benchmarking input format 

Before converting the VCF, please pre-process the output from variant calling pipeline as given [here](https://github.com/broadinstitute/CTAT-benchmarking/tree/master/CTAT-mutation-benchmarking/genome-based-benchmarking#pipeline-output-post-processing). 

1. Run the script
   ```
   python convert_vcf_to_std.py \
   --fname input_variants.vcf \
   --ref ref_chr_GrCh37_high_conf_onlySNPs.vcf \
   --bam alignment.bam \
   --out snps_pipeline.txt \
   --filter_indels
   ```
   Input files : 
   ```
   input_variants.vcf - VCF output file from variant calling pipeline
   ref_chr_GrCh37_high_conf_onlySNPs.vcf - VCF file containing high confidence SNPs
   alignment.bam - BAM file containing alignment between the reads and the reference
   --filter_indels filters the indels and mutiple SNPs
   ```
   
   Output file format :
   ```
   The output file contains the following fields:
   Chr - Chromosome number
   Pos - Position of the SNP
   SNP - Alternative SNP predicted by a prediction algorithm 
   Ref_SNP - SNP given in the reference
   Depth - Coverage of the SNP
   ```
   
