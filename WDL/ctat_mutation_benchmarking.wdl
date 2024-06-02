version 1.0

workflow ctat_mutations_benchmarking {
    input {
        String sample_id
        
        File truth_vcf
        Array[File] predicted_vcf
        File bam
        File? bai

        File? high_conf_regions_bed
        File rna_edit_vcf

        File? exome_bam
        File? exome_bai
        Int? min_exome_depth

        Boolean snvs_only = false
        Boolean indels_only = false
        
        String memory = "18G"
        Int cpu = 8
        Int extra_disk_space = 20

        String docker = "trinityctat/ctat_mutations_benchmarking:latest"
    }

    call benchmark_variants {
        input:
            sample_id=sample_id,
            truth_vcf=truth_vcf,
            predicted_vcf=predicted_vcf,
            bam=bam,
            bai=bai,
            exome_bam=exome_bam,
            exome_bai=exome_bai,
            min_exome_depth=min_exome_depth,
            snvs_only=snvs_only,
            indels_only=indels_only,
            memory=memory,
            extra_disk_space=extra_disk_space,
            cpu=cpu,
            high_conf_regions_bed=high_conf_regions_bed,
            rna_edit_vcf=rna_edit_vcf,
            docker = docker
    }

    output {
        File outputs_tar_gz = benchmark_variants.outputs_tar_gz
    }

    parameter_meta {
        predicted_vcf:{help:"VCF output file from variant calling pipeline"}
        truth_vcf:{help:"VCF file containing high confidence SNPs"}
        bam:{help:"BAM file containing alignment between the reads and the reference"}
    }
}

task benchmark_variants {
    input {

        String sample_id
        
        Array[File] predicted_vcf
        File truth_vcf
        File bam
        File? bai


        File? exome_bam
        File? exome_bai
        Int? min_exome_depth

        Boolean snvs_only = false
        Boolean indels_only = false

        
        File? high_conf_regions_bed
        File rna_edit_vcf

        Int cpu
        Int extra_disk_space
        String memory
        String docker
    }

    String bmark_output_dir = sample_id + "-benchmarking"
    
    command <<<
        set -ex

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

        bmarking_output_dir=~{bmark_output_dir}
        if [[ "~{exome_bam}" != "" ]]; then
            bmarking_output_dir="${bmarking_output_dir}.by_exome"
        fi

        if [[ "~{snvs_only}" == "true" ]]; then
           bmarking_output_dir="${bmarking_output_dir}.snvs_only"
        fi

        if [[ "~{indels_only}" == "true" ]]; then
          bmarking_output_dir="${bmarking_output_dir}.indels_only"
        fi
        
        python /software/ctat-mutations-benchmarking/BENCHMARK_variant_calling.py \
          --pred_vcf ~{sep=' ' predicted_vcf} \
          --truth_vcf ~{truth_vcf} \
          --pred_bam ~{bam} \
          ~{"--exome_bam " + exome_bam} \
          ~{"--min_exome_depth " + min_exome_depth } \
          ~{"--restrict_regions_bed " + high_conf_regions_bed} \
          --rna_editing ~{rna_edit_vcf} \
          --output_dir ~{bmark_output_dir} \
          ~{true='--snvs_only' false='' snvs_only} \
          ~{true='--indels_only' false='' indels_only} \
          --cpu ~{cpu}

         tar -zcf  ~{bmark_output_dir}.tar.gz  ~{bmark_output_dir}
        
    >>>

    runtime {
        disks: "local-disk " + ceil(size(bam, "GB")  + extra_disk_space + size(truth_vcf, "GB")) + " HDD"
        docker: docker
        memory: memory
        preemptible: 2
        cpu : cpu
    }

    output {
        File outputs_tar_gz = select_first(glob("~{sample_id}-benchmarking*tar.gz"))
    }
}
