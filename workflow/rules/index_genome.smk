rule index_genome:
    params:
        {workflow.cores}
    input:
        fna = join(genome_dir, genome_name + '.' + genome_ext)
    output:
        idx = join(genome_dir, genome_name + '.' + genome_ext + ".bwt")
    shell:
        "bwa index {input.fna}"
