rule fastqc:
    params:
        outdir = fastqc
    input:
        join(fq_dir, "{sample}." + fq_ext)
    output:
        join(fastqc, "{sample}_fastqc.html")
    shell:
        "fastqc {input} -o {params.outdir}"
