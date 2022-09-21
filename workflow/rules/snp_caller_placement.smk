rule snp_caller_placement:
    threads:
        workflow.cores
    params:
        perl_script  = "scripts/snp_caller_placement.pl",
        vcf_file_dir = bcftools,
        placed_file  = join(snp_caller_placement, genome_name + "_placed.infile"),
        output_dir   = snp_caller_placement,
    input:
        position_file = join(metadata, genome_name) + '.' + 'positions.txt',
    output:
        #placed_file  = join(snp_caller_placement, genome_name + "_placed.infile"),
        op_file = join(snp_caller_placement, "{sample}.txt")
    shell:
        "perl {params.perl_script} --chrom {genome_name} --chromv 1 --qual 60 --dp 15 --dp4 75 \
        --dpmax 5000 --mq 60 --noindels --noheader --showfiltered --positions {input.position_file} \
        --vcf {wildcards.sample},{params.vcf_file_dir}/{genome_name}_{wildcards.sample}.all.vcf.gz \
        --dir {params.output_dir} --phylip {genome_name}.isolate.infile --verbose 2 --refilter \
        -b 1 --cpus 1 >> {params.placed_file}"
