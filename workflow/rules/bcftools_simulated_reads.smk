rule bcftools_simulated_reads:
    threads:
        5 #workflow.cores
    params:
        tmp_file = join(simulated_read_alignment, genome_name + "_" + "tmp.tmp")
    input:
        ref = join(genome_dir, genome_name) + '.' + genome_ext,
        bam = join(simulated_read_alignment, "{sample}.bam"),
    output:
        vcf = join(bcftools_simulated_reads, genome_name + '_' + "{sample}.all.vcf.gz")
    shell:
        """
            bcftools mpileup -Ou -f {input.ref} {input.bam} --annotate 'FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR' | \
                        bcftools call --threads {workflow.cores} -c | \
                        bcftools filter -sMixed -e '(DP4[0]+DP4[1])>10 & (DP4[2]+DP4[3])>10' - | \
                        bcftools filter -sDepth -e 'FORMAT/DP<10' - | \
                        bcftools filter -sLowQual -g3 -G10 -e '%QUAL<30 || RPB<0.1' - | \
                        bgzip > {output.vcf} && tabix -p vcf {output.vcf} && > {params.tmp_file}
        """

