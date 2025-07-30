version 1.0

struct RuntimeEnvironment {
    String docker_image
    Int preemptible
    Int boot_disk_size
    Int disk_space
    Int cpu
    Int mem 
}

workflow CtatLRMutationWorkflow {
    meta {
        version: "v0.1"
        author: "Qian Qin"
    }

    input {
        File? fastq
        File? ubam

        File? SplitNCigarBam
        File? SplitNCigarBai
        File? AlignedBam
        File? AlignedBai

        String sample_id
        File fasta
        File fasta_idx
        File fasta_dict

        File gtf
        String? readType = "PacBioIsoSeqUb"

        Int cpu = 24
        Int scatterCount = 24
        Int scatter_mem = 36
        Int scatter_cpu = 12
        Int mem = 50
        Int disk_space
        Int preemptible
        Int boot_disk_size

        Boolean is_split_chrom_for_deepvariant = true
        Boolean run_clair3_rna = true
        Boolean run_longcallR = true
        Boolean run_clair3_rna_mix = true
        Boolean run_clair3_mix = true
        Boolean run_deepvariant = true

        File monitoring_script = "gs://broad-dsde-methods-monitoring/cromwell_monitoring_script.sh"
    }

    RuntimeEnvironment runtime_environment = {
        'docker_image': "trinityctat/longgf", 'preemptible': preemptible, 'boot_disk_size': boot_disk_size, 'disk_space': disk_space, 
        'cpu': cpu, 'mem': mem
    }

    RuntimeEnvironment runtime_environment2 = {
        'docker_image': "trinityctat/flagcorrection", 'preemptible': preemptible, 'boot_disk_size': boot_disk_size, 'disk_space': disk_space, 
        'cpu': cpu, 'mem': mem
    }

    RuntimeEnvironment runtime_environment4 = {
        'docker_image': "hkubal/clair3-rna:latest", 'preemptible': preemptible, 'boot_disk_size': boot_disk_size, 'disk_space': disk_space, 
        'cpu': cpu, 'mem': mem
    }

    RuntimeEnvironment runtime_environment5 = {
        'docker_image': "hkubal/clair3:latest", 'preemptible': preemptible, 'boot_disk_size': boot_disk_size, 'disk_space': disk_space, 
        'cpu': cpu, 'mem': mem
    }

    RuntimeEnvironment runtime_environment6 = {
        'docker_image': "qianqin/longcallr", 'preemptible': preemptible, 'boot_disk_size': boot_disk_size, 'disk_space': disk_space, 
        'cpu': cpu, 'mem': mem
    }

    String scripts_path = "/usr/local/src/ctat-mutations/src"
    String docker = "trinityctat/ctat_mutations:4.3.0"

    if ( defined(ubam) && (!defined(SplitNCigarBam)) && (!defined(AlignedBam)) && (!defined(fastq)) ) {
        call samtools_fastq {
            input:
                sample_id = sample_id,
                ubam = ubam,
                runtime_environment = runtime_environment,
                monitoring_script = monitoring_script
        }
    }

	if ((!defined(SplitNCigarBam)) && (!defined(SplitNCigarBai)) && (!defined(AlignedBam)) && (!defined(AlignedBai)) ) {
	    File input_fastq = select_first([fastq, samtools_fastq.fastq])
    }

    if ( !defined(input_fastq) && (!defined(SplitNCigarBam)) && (!defined(AlignedBam)) ) {
        call raise_exception as error_input_data  { 
            input:
                msg = "No FASTQ/SplitNCigarBam/AlignedBam input",
                runtime_environment = runtime_environment
        }
    }

    if ( defined(input_fastq) ) {
        call minimap2Task {
            input:
                 fastq = input_fastq,
                 sample_id = sample_id,
                 gtf   = gtf,
                 fasta = fasta,
                 readType = readType,
                 runtime_environment = runtime_environment,
                 monitoring_script = monitoring_script
        }
    }

    File input_bam = select_first([minimap2Task.bam, AlignedBam, SplitNCigarBam])
    File input_bai = select_first([minimap2Task.bai, AlignedBai, SplitNCigarBai])
    
    if ((!defined(SplitNCigarBam)) && (!defined(SplitNCigarBai))) {    
      if ( run_clair3_rna ) {
          call clairRNATask {
              input:
                  input_bam=input_bam,
                  input_bai=input_bai,
                  sample_id=sample_id,
                  assembly=fasta,
                  fai=fasta_idx,
                  runtime_environment=runtime_environment4
          }
      }

      if ( run_longcallR ) {
          call longcallrTask {
              input:
                  input_bam=input_bam,
                  input_bai=input_bai,
                  sample_id=sample_id,
                  assembly=fasta,
                  fai=fasta_idx,
                  runtime_environment=runtime_environment6
          }
      }      
    }

    if ( run_clair3_rna_mix || run_clair3_mix || run_deepvariant ) {
        if ( (!defined(SplitNCigarBam)) && (!defined(SplitNCigarBai)) ) {
            call SplitNCigarLongReads {
                input:
                    input_bam = input_bam,
                    input_bam_index = input_bai,
                    scripts_path = scripts_path,
                    docker = docker,
                    preemptible = preemptible,
                    monitoring_script = monitoring_script
            }
        }
        File split_ncigar_bam = select_first([SplitNCigarLongReads.bam, input_bam])
        File split_ncigar_bai = select_first([SplitNCigarLongReads.bai, input_bai])

        if ( run_clair3_rna_mix ) {
            call clairRNATask as clair_rna_mix {
                input:
                    input_bam=split_ncigar_bam,
                    input_bai=split_ncigar_bai,
                    sample_id=sample_id,
                    assembly=fasta,
                    fai=fasta_idx,
                    runtime_environment=runtime_environment4
            }
        }

        if ( run_deepvariant ) {
           RuntimeEnvironment runtime_environment3 = {
               'docker_image': "google/deepvariant", 'preemptible': preemptible, 'boot_disk_size': boot_disk_size, 'disk_space': disk_space, 
               'cpu': scatter_cpu, 'mem': scatter_mem
           }

           if (is_split_chrom_for_deepvariant) {
               Array[String] chroms = ["chr1", "chr2", "chr3", "chr4", "chr5",
                                       "chr6", "chr7", "chr8", "chr9", "chr10",
                                       "chr11", "chr12", "chr13", "chr14", "chr15",
                                       "chr16", "chr17", "chr18", "chr19", "chr20",
                                       "chr21", "chr22", "chrX"]
	       scatter (interval in chroms) {
                   call DeepVariant {
                       input:
                           bam = split_ncigar_bam,
                           bai = split_ncigar_bai,
                           sample_id = sample_id,
                           fasta = fasta,
                           fasta_idx = fasta_idx,
                           region = interval,
                           runtime_environment = runtime_environment3,
                           monitoring_script = monitoring_script
                   }
	           File DeepVariantVcf = DeepVariant.vcf
	           File DeepVariantVcfIdx = DeepVariant.tbi
	       }

               call MergeVCFs {
                   input:
                       input_vcfs = DeepVariantVcf,
                       input_vcfs_indexes = DeepVariantVcfIdx,
                       sample_id = sample_id
               }
           }

           if (!is_split_chrom_for_deepvariant) {
               call DeepVariant as deepvariant_all_chrom {
                   input:
                       bam = split_ncigar_bam,
                       bai = split_ncigar_bai,
                       sample_id = sample_id,
                       fasta = fasta,
                       fasta_idx = fasta_idx,
                       runtime_environment = runtime_environment3,
                       monitoring_script = monitoring_script
               }
           }
        }

        if ( run_clair3_mix ) {
           call Clair3Mix {
               input:
                   input_bam= split_ncigar_bam,
                   input_bai= split_ncigar_bai,
                   sample_id = sample_id,
                   assembly=fasta,
                   fai=fasta_idx,
                   runtime_environment=runtime_environment5,
                   monitoring_script = monitoring_script
           }
        }
    }
    
    output {
        File? minimap2_bam = minimap2Task.bam
        File? minimap2_bai = minimap2Task.bai
        File? minimap2_log = minimap2Task.monitoring_log
        File? split_bam = split_ncigar_bam
        File? split_bai = split_ncigar_bai

        File? clair_rna_vcf = clairRNATask.vcf
        File? clair_rna_mix_vcf = clair_rna_mix.vcf

        Array[File]? clair_dna_mix_vcfdir = Clair3Mix.outdir
        File? clair_dna_mix_log = Clair3Mix.log

        File? deepvariant_merged_vcf = select_first([MergeVCFs.vcf, deepvariant_all_chrom.vcf])
        File? deepvariant_merged_tbi = select_first([MergeVCFs.tbi, deepvariant_all_chrom.tbi])

        File? longcallr_vcf = longcallrTask.vcf
        File? longcallr_bam = longcallrTask.bam
    }
}

task samtools_fastq {
    input {
        String sample_id
        File? ubam
        RuntimeEnvironment runtime_environment

        File monitoring_script 
    }

    command <<< 
        bash ~{monitoring_script} > ~{sample_id}_monitoring.log &

        samtools fastq -@ ~{runtime_environment.cpu} ~{ubam} | gzip - > ~{sample_id}.fastq.gz
    >>>

    output {
        File fastq = "~{sample_id}.fastq.gz"
        File monitoring_log = "~{sample_id}_monitoring.log"
    }

    runtime {
        disks: "local-disk ~{runtime_environment.disk_space} HDD"
        memory: "~{runtime_environment.mem} GB"
        cpu: runtime_environment.cpu
        preemptible: runtime_environment.preemptible
        bootDiskSizeGb: runtime_environment.boot_disk_size
        docker: runtime_environment.docker_image
    }
}

task longcallrTask {
    input {
        File input_bam
        File input_bai
        String sample_id
        File assembly
        File fai

        RuntimeEnvironment runtime_environment
    }

    command <<<
        longcallR --bam-path ~{input_bam} --ref-path ~{assembly} --platform hifi --preset hifi-masseq --output ~{sample_id}
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
        File vcf = "~{sample_id}.vcf"
        File bam = "~{sample_id}.phased.bam"
    }
}

task minimap2Task {
    input {
        File? fastq
        String sample_id
        String? readType
        File fasta
        File gtf

        RuntimeEnvironment runtime_environment
        File monitoring_script 
    }

    command <<<
        bash ~{monitoring_script} > ~{sample_id}_monitoring.log &

        # copy from https://github.com/broadinstitute/MDL-workflows/blob/main/LR-tools/minimap2_LR/minimap2_LR.wdl
        if [ "~{readType}" == "SplicedLongReads" ]; then
            minimap2_preset="splice"
        elif [ "~{readType}" == "ONTDirectRNA" ]; then
            minimap2_preset="splice -uf -k14"
        elif [ "~{readType}" == "PacBioIsoSeq" ]; then
            minimap2_preset="splice:hq "
        elif [ "~{readType}" == "PacBioIsoSeqUf" ]; then
            minimap2_preset="splice:hq -uf"
        elif [ "~{readType}" == "PacBioIsoSeqUb" ]; then
            minimap2_preset="splice:hq -ub"
        else
            echo "Invalid readType: ~{readType}"
            exit 1
        fi

        paftools.js gff2bed ~{gtf} > junc.bed
        minimap2 --secondary=no -C5 -t ~{runtime_environment.cpu} -ax ${minimap2_preset} --junc-bed junc.bed ~{fasta} ~{fastq} | samtools view -bSh -F 2308 - > ~{sample_id}.bam
        samtools sort -@ ~{runtime_environment.cpu} -o ~{sample_id}_sort.bam ~{sample_id}.bam
        samtools index -@ ~{runtime_environment.cpu} ~{sample_id}_sort.bam
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
        File bam = "~{sample_id}_sort.bam"
        File bai = "~{sample_id}_sort.bam.bai"
        File monitoring_log = "~{sample_id}_monitoring.log"
    }
}

task SplitNCigarLongReads {
    # use Brian's script instead
    input {
        File input_bam
        File input_bam_index
        String scripts_path
        
        String docker
        Int preemptible
        File monitoring_script
        
    }

    String output_bam_filename = basename(input_bam, ".bam") + ".splitNcigar.bam"
    
    command <<<
        set -ex

        ~{scripts_path}/cigar_N_splitter.py ~{input_bam} split_N.bam

        samtools sort split_N.bam -o ~{output_bam_filename}

        samtools index  ~{output_bam_filename}
        
    >>>


    output {
        File bam = "~{output_bam_filename}"
        File bai = "~{output_bam_filename}.bai"
    }
        
    runtime {
        disks: "local-disk " + ceil((size(input_bam, "GB") + 10) * 10 ) + " SSD"
        docker: docker
        memory: "8GB"
        preemptible: preemptible
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
        File monitoring_script
    }

    command <<<
        ulimit -u 10000 
        bash ~{monitoring_script} > ~{sample_id}_monitoring.log &

        ##grep -v "@" ~{region} | cut -f 1-3 | awk -v OFS="\t" '{print $1,$2-1,$3+1}' > region.bed

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
        File log = "~{sample_id}_monitoring.log"
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

task Clair3Mix {
    input {
        File input_bam
        File input_bai
        String sample_id
        File assembly
        File fai
        RuntimeEnvironment runtime_environment
        File monitoring_script
    }

    command <<<
        bash ~{monitoring_script} > ~{sample_id}_monitoring.log &
        /opt/bin/run_clair3.sh \
          --bam_fn=~{input_bam} \
          --ref_fn=~{assembly} \
          --threads=~{runtime_environment.cpu} \
          --platform="hifi" \
          --model_path=/opt/models/hifi_revio \
          --output=~{sample_id}_clair3
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
        Array[File]? outdir = glob("~{sample_id}_clair3/*")
        File log = "~{sample_id}_monitoring.log"
    }
}

task clairRNATask {
    input {
        File? input_bam
        File? input_bai
        String sample_id
        File assembly
        File fai

        RuntimeEnvironment runtime_environment
    }

    command <<<
        /opt/bin/run_clair3_rna -B ~{input_bam} -R ~{assembly} -o ~{sample_id} -t ~{runtime_environment.cpu} --min_mq 30  -p hifi_mas        
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
        File vcf = "~{sample_id}/output.vcf.gz"
    }
}

task raise_exception {
    input {
        String msg
        RuntimeEnvironment runtime_environment
    }
    command {
        echo -e "\n* Error: ${msg}\n" >&2
        exit 2
    }
    output {
        String error_msg = '${msg}'
    }
    runtime {
        maxRetries : 0
        cpu : 1
        memory : '2 GB'
        time : 4
        disks : 'local-disk 10 SSD'
        docker : runtime_environment.docker_image
    }
}