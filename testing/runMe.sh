../BENCHMARK_variant_calling.py \
    --pred_vcf chr18.PRED.raw.vcf chr18.PRED.init.vcf chr18.PRED.annot_n_filtered.vcf \
    --truth_vcf chr18.REF.vcf \
    --pred_bam chr18.bam \
    --restrict_regions_bed chr18.restrict_regions.bed \
    --dbsnp chr18.dbsnp.vcf \
    --rnaediting chr18.rnaediting.vcf \
    --remove_indels \
    --remove_intersect rnaediting \
    --output_dir "bmark.outdir"


