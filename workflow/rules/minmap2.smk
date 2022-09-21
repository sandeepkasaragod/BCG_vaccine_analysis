rule minmap:
    threads:
        workflow.cores
    params:
        tempf="temp.{sample}"
    input:
        ref = join(genome_dir, genome_name) + '.' + genome_ext,
        reads = join(fq_dir, "{sample}.fastq.gz"),
    output:
        bam = join(minmap_alignment, "alignment", "{sample}_" + genome_name + ".bam"),
        alignment_summary = join(minmap_alignment, "{sample}_summary.txt")
    shell:
        "minimap2 -t {threads} -ax map-ont {input.ref} {input.reads} | samtools view -bS - | samtools sort -T {params.tempf} - | samtools rmdup - - > {output.bam} \
                && samtools index {output.bam} && samtools flagstat {output.bam} > {output.alignment_summary}"
