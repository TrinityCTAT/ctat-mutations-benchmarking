


snvs_only:
	../BENCHMARK_variant_calling.py \
    --pred_vcf chr18.PRED.raw.vcf chr18.PRED.init.vcf chr18.PRED.annot_n_filtered.vcf \
    --truth_vcf chr18.REF.vcf \
    --pred_bam chr18.bam \
    --restrict_regions_bed chr18.restrict_regions.bed \
    --rna_editing chr18.rnaediting.vcf.gz \
    --snvs_only \
    --output_dir "bmark.outdir.snvs_only"


indels_only:
	../BENCHMARK_variant_calling.py \
    --pred_vcf chr18.PRED.raw.vcf  \
    --truth_vcf chr18.REF.vcf \
    --pred_bam chr18.bam \
    --restrict_regions_bed chr18.restrict_regions.bed \
    --rna_editing chr18.rnaediting.vcf.gz \
    --indels_only \
    --output_dir "bmark.outdir.indels_only"



clean:
	./cleanMe.sh

