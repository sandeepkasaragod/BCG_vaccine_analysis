rule iqtree:
        threads:
                workflow.cores
        params:
                script_path     = "scripts/placed_file_to_fasta.py",
                output_dir      = iqtree,
                infile          = join(snp_caller_placement, genome_name),
                placement_fasta = join(iqtree, genome_name + '_placed.infile.fasta')
        input:
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
