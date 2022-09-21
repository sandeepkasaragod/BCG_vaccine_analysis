rule wgsim:
    threads:
        1 #workflow.cores
    input:
        fna = join(config['reference_strains'], "{sample}.fna"),
    output:
        r1 = join(wgsim, "{sample}_1.fastq"),
        r2 =  join(wgsim, "{sample}_2.fastq")
    conda:
        "envs/wgsim.yml"
    shell:
        """
                wgsim {input.fna} {output.r1} {output.r2}
        """
