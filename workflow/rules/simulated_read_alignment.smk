rule simulated_read_alignment:
    threads:
        workflow.cores
    params:
        tempf='temp.{sample}'
    input:
        ref = join(genome_dir, genome_name) + '.' + genome_ext,
        r1 = join(wgsim, "{sample}_1.fastq"),
        r2 = join(wgsim, "{sample}_2.fastq"),
    output:
        bam = join(simulated_read_alignment, "{sample}.bam"),
        alignment_summary = join(simulated_read_alignment, "{sample}_summary.txt"),
    shell:
        "bwa mem -t {threads} {input.ref} {input.r1} {input.r2} | samtools view -bS - | samtools sort -T {params.tempf} - | samtools rmdup - - > {output.bam} \
         && samtools index {output.bam} && samtools flagstat {output.bam} > {output.alignment_summary}"
