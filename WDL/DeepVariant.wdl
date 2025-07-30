version 1.0

struct RuntimeEnvironment {
    String docker_image
    Int preemptible
    Int boot_disk_size
    Int disk_space
    Int cpu
    Int mem 
}


workflow DeepVariant_workflow {

  input {
  
    String sample_id
    File ref_genome_fasta
    File ref_genome_fasta_idx
    
    File AlignedBam
    File AlignedBai
    
    Int cpu = 24
    Int scatterCount = 24
    Int scatter_mem = 36
    Int scatter_cpu = 12
    Int mem = 50
    Int disk_space = 200
    Int preemptible = 0
    Int boot_disk_size = 200

    Boolean scatter_chromosomes = false
  }
  
  RuntimeEnvironment runtime_env = {
      'docker_image': "google/deepvariant",
      'preemptible': preemptible,
      'boot_disk_size': boot_disk_size,
      'disk_space': disk_space, 
      'cpu': scatter_cpu,
      'mem': scatter_mem
    }

    if (scatter_chromosomes) {
      Array[String] chroms = ["chr1", "chr2", "chr3", "chr4", "chr5",
                              "chr6", "chr7", "chr8", "chr9", "chr10",
                              "chr11", "chr12", "chr13", "chr14", "chr15",
                              "chr16", "chr17", "chr18", "chr19", "chr20",
                              "chr21", "chr22", "chrX"]
      
      scatter (interval in chroms) {
        call DeepVariant as DVscattered {
          input:
          bam = AlignedBam,
          bai = AlignedBai,
          sample_id = sample_id,
          fasta = ref_genome_fasta,
          fasta_idx = ref_genome_fasta_idx,
          region = interval,
          runtime_environment = runtime_env
        }
	  }

      call MergeVCFs {
        input:
        input_vcfs = DVscattered.vcf,
        input_vcfs_indexes = DVscattered.tbi,
        sample_id = sample_id
      }
    }
        
    if (! scatter_chromosomes) {
      call DeepVariant as deepvariant_all_chrom {
        input:
        bam = AlignedBam,
        bai = AlignedBai,
        sample_id = sample_id,
        fasta = ref_genome_fasta,
        fasta_idx = ref_genome_fasta_idx,
        runtime_environment = runtime_env
      }
    }
  
  
    output {
      File deepvariant_merged_vcf = select_first([MergeVCFs.vcf, deepvariant_all_chrom.vcf])
      File deepvariant_merged_tbi = select_first([MergeVCFs.tbi, deepvariant_all_chrom.tbi])
    }
}

task DeepVariant {
    input {
        File bam
        File bai
        File fasta
        File fasta_idx
        String sample_id
        String? region
        RuntimeEnvironment runtime_environment
    }

    command <<<
        ulimit -u 10000 
    
        python /opt/deepvariant/bin/run_deepvariant.py \
        --model_type PACBIO \
        --ref ~{fasta} \
        --reads ~{bam} \
        --output_vcf ~{sample_id}.vcf.gz \
        --num_shards ~{runtime_environment.cpu} \
        ~{"--regions " + region}
    >>>

    runtime {
        disks: "local-disk ~{runtime_environment.disk_space} HDD"
        memory: "~{runtime_environment.mem} GB"
        cpu: runtime_environment.cpu
        preemptible: runtime_environment.preemptible
        bootDiskSizeGb: runtime_environment.boot_disk_size
        docker: runtime_environment.docker_image
    }

     output {
        File vcf = "~{sample_id}.vcf.gz"
        File tbi = "~{sample_id}.vcf.gz.tbi"
    }
}


task MergeVCFs {
  input {
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String sample_id

    String gatk_path = "/gatk/gatk"
    String docker = "broadinstitute/gatk:latest"
  }

  # Using MergeVcfs instead of GatherVcfs so we can create indices
  # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
  command <<<
    ~{gatk_path} \
    MergeVcfs \
    --INPUT ~{sep=' --INPUT ' input_vcfs} \
    --OUTPUT ~{sample_id}_deepvariant.vcf.gz
  >>>

  output {
    File vcf = "~{sample_id}_deepvariant.vcf.gz"
    File tbi = "~{sample_id}_deepvariant.vcf.gz.tbi"
  }

  runtime {
    memory: "8 GB"
    disks: "local-disk " + 200 + " HDD"
    docker: docker
    preemptible: 3
  }
}