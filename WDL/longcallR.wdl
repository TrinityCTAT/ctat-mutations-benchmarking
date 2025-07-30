version 1.0


workflow longcallrWorkflow {
    input {
        File bam
        File bai
        String sample_id
        File assembly
        File fai
        String docker_image="qianqin/longcallr"
        Int preemptible=3
        Int boot_disk_size=20
        Int disk_space=300
        Int cpu = 10
        Int mem = 40
    }

    call longcallrTask {
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
            mem=mem,
    }

    output {
        File phased_bam = longcallrTask.output_bam
        File phased_vcf = longcallrTask.output_vcf
    }
}


task longcallrTask {
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
        longcallR --bam-path ~{input_bam} --ref-path ~{assembly} --platform hifi --preset hifi-masseq --output ~{sample_id}
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
        File output_vcf = "~{sample_id}.vcf"
        File output_bam = "~{sample_id}.phased.bam"
    }
}