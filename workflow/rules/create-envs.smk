my_envs = ['bowtie2.yaml', 'cutadapt.yaml', 'fastqc.yaml', 'fastqscreen.yaml', 'kraken2.yaml', 'samtools.yaml', 'seqtk.yaml']

rule make_all_envs:
    input:
        expand("created-{name}", name=my_envs)


for env_file in my_envs:
    rule:
        output:
            temp("created-%s" % env_file)
        conda:
            "envs/%s" % env_file
        shell:
            "touch {output}"
