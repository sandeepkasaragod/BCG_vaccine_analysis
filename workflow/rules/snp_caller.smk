rule snp_caller:
    threads:
        workflow.cores
    params:
        outdir       = join(snp_caller),
        bcftools_res = join(bcftools_simulated_reads),
    input:
        meta_data     = join(metadata, "file_list_reference.txt"),
        bcftools_resl = expand(join(bcftools_simulated_reads, genome_name + '_' + "{sample}.all.vcf.gz"), sample=reference_strains),
    output:
        wgsim_pos     = join(snp_caller, "wgsim.positions.txt"),
        fasta_file    = join(snp_caller, "wgsim"),
        position_file = join(metadata, genome_name) + '.' + 'positions.txt',
    shell:
        "perl scripts/snp_caller.pl --chrom {genome_name} --chromv 1 --qual 30 --dp 4 --dp4 75 --dpmax 5000 --af 0 --mq 60 --noindels --noheader \
        --vcfdir {params.bcftools_res} \
                --include {input.meta_data} \
        --dir {params.outdir} \
        --phylip wgsim --verbose 1 --refilter -b 1 --cpus {threads} \
        --include {input.meta_data} && cp {output.wgsim_pos} {output.position_file}"
