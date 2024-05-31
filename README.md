# RNA-seq Variant Calling Benchmarking

## Benchmarking using a truth-set VCF


```

    ../BENCHMARK_variant_calling.py \
        --pred_vcf methodA.vcf [methodB.vcf methodC.vcf ...] \
        --truth_vcf truth.vcf \
        --pred_bam rnaseq_alignments.bam \
        --rna_editing rnaediting.vcf.gz \
        --output_dir "bmark.outdir"

     and optionally:
         [ --restrict_regions_bed trusted_regions.bed ] \
         [ --snvs_only or --indels_only ]

```


>Separately evaluate SNVs or InDels using --snvs_only or --indels_only, respectively.

> the rnaediting.vcf.gz is based on rediportal and provided in the ctat genome lib.

## Benchmarking using a matched whole exome based vcf

First, generate an exome-based VCF using the standard GATK exome variant calling pipeline.

Then, provide the exome vcf as one of the inputs like so:

```

 ../BENCHMARK_variant_calling.py \
    --pred_vcf chr18.PRED.raw.vcf chr18.PRED.init.vcf chr18.PRED.annot_n_filtered.vcf \
    --truth_vcf exome.vcf \
    --pred_bam  rnaseq_alignments.bam \
    --exome_bam exome_dna_alignments.bam \
    --min_exome_depth 10 \
    --rna_editing rnaediting.vcf.gz \
    --output_dir "bmark.outdir.exome_based"

  and optionally:
    [ --snvs_only or --indels_only]

```

