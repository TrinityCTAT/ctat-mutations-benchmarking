version 1.0


workflow ClairWorkflow {
    input {
        File bam
        File bai
        String sample_id
        File assembly
        File fai
        String docker_image="hkubal/clair3-rna:latest"
        Int preemptible=3
        Int boot_disk_size=20
        Int disk_space=300
        Int cpu = 10
        Int mem = 40
    }

    call clairTask {
        input:
            input_bam=bam,
            input_bai=bai,
            sample_id=sample_id,
            docker_image=docker_image,
            assembly=assembly,
            fai=fai,
            preemptible=preemptible,
            boot_disk_size=boot_disk_size,
            disk_space=disk_space,
            cpu=cpu,
            mem=mem
    }

    output {
        File phased_vcf = clairTask.output_vcf_gz
    }
}


task clairTask {
    input {
        File input_bam
        File input_bai
        String sample_id
        File assembly
        File fai
        String docker_image
        Int preemptible=3
        Int boot_disk_size=10
        Int disk_space=20
        Int cpu = 10
        Int mem = 64
    }

    command <<<
        /opt/bin/run_clair3_rna -B ~{input_bam} -R ~{assembly} -o ~{sample_id} -t ~{cpu} --min_mq 30  -p hifi_mas        
    >>>

    runtime {
        disks: "local-disk ~{disk_space} HDD"
        memory: "~{mem} GB"
        cpu: cpu
        preemptible: preemptible
        bootDiskSizeGb: boot_disk_size
        docker: docker_image
    }

    output {
        File output_vcf_gz = "~{sample_id}/output.vcf.gz"
    }
}