rule multiqc:
    params:
        qcdir = fastqc,
        outdir = multiqc
    input:
        expand(join(fastqc, "{sample}_fastqc.html"), sample=fq_files)
    output:
        outfile = join(multiqc, "multiqc_report.html"),
    shell:
        "multiqc {params.qcdir} -o {params.outdir}"
