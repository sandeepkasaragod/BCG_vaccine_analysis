import sys
from os.path import join as join
import pandas as pd
import os

configfile: 'config.yml'
metadata = "metadata"

#creating mandate dirs
os.makedirs(config['output_dir'], exist_ok=True)
os.makedirs(join(config['output_dir'], 'bcftools_simulated_reads'), exist_ok=True)

#result dir paths
fastqc = join(config['output_dir'], "fastqc") 
multiqc = join(config['output_dir'], "multiQC")

minmap_alignment = join(config['output_dir'], "minimap2")

bcftools = join(config['output_dir'], "bcftools") 

wgsim = join(config['output_dir'], "simulated_reads") 

simulated_read_alignment = join(config['output_dir'], "simulated_read_alignment")
bcftools_simulated_reads = join(config['output_dir'], "bcftools_simulated_reads")

snp_caller = join(config['output_dir'], "snp_caller")
snp_caller_placement = join(config['output_dir'], "snp_caller_placement")

iqtree = join(config['output_dir'], "iqtree")

#load input fastq files
fq_dir = config['fq_dir']
fq_files = list(set([x.split('.')[0] for x in os.listdir(fq_dir)]))
fq_ext = config['file_extension']


#loading reference straisn
reference_strains = list(set([x.split('.')[0] for x in os.listdir(config['reference_strains'])]))

#load genome file
genome_dir = config['ref_genome_dir']
genome_name = config['ref_genome_name']
genome_ext = config['ref_genome_ext']
genome_version = config['ref_genome_version']

#main rule
rule all:
	input:
		join(genome_dir, genome_name + '.' + genome_ext + ".bwt"), #bwa_index
		expand(join(config['fq_dir'], "{sample}." + fq_ext), sample=fq_files), #fastqc
		join(multiqc, "multiqc_report.html"), #multiqc
		expand(join(minmap_alignment, "alignment", "{sample}_" + genome_name + ".bam"), sample=fq_files), #minmap2
		expand(join(bcftools, genome_name + "_" + "{sample}.all.vcf.gz"), sample=fq_files), #bcftools
		expand(join(wgsim, "{sample}_1.fastq"), sample=reference_strains), #generating simulated reads wgsim
		expand(join(simulated_read_alignment, "{sample}.bam"), sample=reference_strains), #align the simulated reads
		expand(join(bcftools_simulated_reads, genome_name + '_' + "{sample}.all.vcf.gz"), sample=reference_strains), #bcftools for simulated reads
		join(snp_caller, "wgsim.positions.txt"), #snpcaller for reference reads
		expand(join(snp_caller_placement, "{sample}.txt"), sample=fq_files),
		join(iqtree, genome_name + '.treefile'), #iqtree


#indexing the genome
rule index_genome:
	params:
		{workflow.cores}
	input:
		fna = join(genome_dir, genome_name + '.' + genome_ext)
	output:
		idx = join(genome_dir, genome_name + '.' + genome_ext + ".bwt")
	shell:
		"bwa index {input.fna}"


#qc using fastqc	
rule fastqc:
        params:
                outdir = fastqc
        input:
                join(fq_dir, "{sample}." + fq_ext)
        output:
                join(fastqc, "{sample}_fastqc.html")
        shell:
                "fastqc {input} -o {params.outdir}"


#combining qc using multiqc
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


#alignment using Minmap2
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


#bcftools 
rule bcftools:
    threads:
        1 #workflow.cores
    input:
        ref = join(genome_dir, genome_name) + '.' + genome_ext,
        bam = join(minmap_alignment, "alignment", "{sample}_" + genome_name + ".bam"),
    output:
        vcf = join(bcftools, genome_name + "_" + "{sample}.all.vcf.gz")
    shell:
        """
            bcftools mpileup -Ou -f {input.ref} {input.bam} --annotate 'FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR' | \
                        bcftools call --threads 1 -c | \
                        bcftools filter -sMixed -e '(DP4[0]+DP4[1])>10 & (DP4[2]+DP4[3])>10' - | \
                        bcftools filter -sDepth -e 'FORMAT/DP<10' - | \
                        bcftools filter -sLowQual -g3 -G10 -e '%QUAL<30 || RPB<0.1' - | \
                        bgzip > {output.vcf} && tabix -p vcf {output.vcf}
        """
#-cv to show only variants called


#generating simulated reads
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


#aligning the simulated reads using bwa alignment
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


#bcftool calling on simulated reads
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

#site calling for reference reads
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

#--phylip wgsim --verbose 1 --refilter -b 1 --cpus {threads}"


#site calling for samples
rule snp_caller_placement:
	threads:
		workflow.cores
	params:
		perl_script  = "scripts/snp_caller_placement.pl",
		vcf_file_dir = bcftools,
		placed_file  = join(snp_caller_placement, genome_name + "_placed.infile"),
		output_dir   = snp_caller_placement,
		
	input:
		#each_sample  = expand(join(bcftools, genome_name + "_" + "{sample}.all.vcf.gz"), sample=fq_files),
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

#tree construction
rule iqtree:
        threads:
                workflow.cores
        params:
                script_path     = "scripts/placed_file_to_fasta.py",
                output_dir      = iqtree,
                infile          = join(snp_caller_placement, genome_name),
                placement_fasta = join(iqtree, genome_name + '_placed.infile.fasta')
        input:
                #infile     = join(snp_caller_placement, genome_name), # this portion is added in the script + '_placed.infile'),
                snp_caller  = join(snp_caller, 'wgsim'),
        output:
                fasta   = join(iqtree, genome_name + '_placed.infile.fasta'),
                outtree = join(iqtree, genome_name + '.treefile'),
                info    = join(iqtree, genome_name + '.log'),
                ckp     = join(iqtree, genome_name  + '.ckp.gz'),
        log:
                'logs/run_iqtree.log'
        shell:
                """
                        python {params.script_path} -i {params.infile} -o {params.output_dir} -s {input.snp_caller} -g {genome_name} -v {genome_version} && iqtree -T {threads} -mset raxml -mfreq F -B 1000 --prefix {params.output_dir} -s {params.placement_fasta} > {log} 2>&1
                """

