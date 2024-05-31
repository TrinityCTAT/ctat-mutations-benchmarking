version 1.0

workflow ctat_mutations_benchmarking {
    input {
        File truth_vcf
        Array[File] predicted_vcf
        File bam
        String memory = "18G"
        Int cpu = 8
        Int extra_disk_space = 20
        File high_conf_regions_bed
        File rna_edit_vcf
        String docker = "trinityctat/ctat_mutation_benchmarking:latest"
    }

    call benchmark_variants {
        input:
            truth_vcf=truth_vcf,
            predicted_vcf=predicted_vcf,
            bam=bam,
            memory=memory,
            extra_disk_space=extra_disk_space,
            cpu=cpu,
            high_conf_regions_bed=high_conf_regions_bed,
            rna_edit_vcf=rna_edit_vcf,
            docker = docker
    }

    output {
        Array[File] outputs = benchmark_variants.outputs
    }

    parameter_meta {
        predicted_vcf:{help:"VCF output file from variant calling pipeline"}
        truth_vcf:{help:"VCF file containing high confidence SNPs"}
        bam:{help:"BAM file containing alignment between the reads and the reference"}
    }
}

task benchmark_variants {
    input {
        Array[File] predicted_vcf
        File truth_vcf
        File bam
        File high_conf_regions_bed
        File rna_edit_vcf
        Int cpu
        Int extra_disk_space
        String memory
        String docker
    }

    command <<<
        set -e

        mkdir benchmarking

        python <<CODE
        import os
        import subprocess
        if not os.path.exists('~{bam}.bai'):
            subprocess.check_call(['samtools', 'index', '~{bam}'])
        vcfs = []
        vcfs.append('~{rna_edit_vcf}')
        for vcf in vcfs:
            if not os.path.exists(vcf + '.csi'):
                subprocess.check_call(['bcftools', 'index', vcf])
        CODE

        python /software/CTAT-benchmarking/CTAT-mutation-benchmarking/BENCHMARK_variant_calling.py \
        --pred_vcf ~{sep=' ' predicted_vcf} \
        --truth_vcf ~{truth_vcf} \
        --pred_bam ~{bam} \
        --restrict_regions ~{high_conf_regions_bed} \
        --rnaediting ~{rna_edit_vcf} \
        --remove_intersect rnaediting \
        --output_dir benchmarking \
        --remove_indels \
        --cpu ~{cpu}
    >>>

    runtime {
        disks: "local-disk " + ceil(size(bam, "GB")  + extra_disk_space + size(truth_vcf, "GB")) + " HDD"
        docker: docker
        memory: memory
        preemptible: 2
        cpu : cpu
    }

    output {
        Array[File] outputs = glob('benchmarking/*')
    }
}
