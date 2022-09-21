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


include: join("workflow", "rules", "index_genome.smk")
include: join("workflow", "rules", "fastqc.smk")
include: join("workflow", "rules", "multiqc.smk")
include: join("workflow", "rules", "minmap2.smk")
include: join("workflow", "rules", "bcftools.smk")
include: join("workflow", "rules", "wgsim.smk")
include: join("workflow", "rules", "simulated_read_alignment.smk")
include: join("workflow", "rules", "bcftools_simulated_reads.smk")
include: join("workflow", "rules", "snp_caller.smk")
include: join("workflow", "rules", "snp_caller_placement.smk")
include: join("workflow", "rules", "iqtree.smk")
