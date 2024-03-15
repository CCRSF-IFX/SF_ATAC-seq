
rule fastqc:
    input: R1 = rawdir + "/" + "{sample}_R1_001.fastq", R2 = rawdir + "/" + "{sample}_R2_001.fastq"
    output: "Sample_{sample}/{sample}_R1_001_fastqc.html"
    params: prefix = "Sample_{sample}/"
    log: "logs/{sample}_fastqc.log"
    conda: "envs/fastqc.yaml"
    threads: 16
    resources: mem_mb=config['medium']['mem'], time=config['medium']['time'], partition=config['medium']['partition']
    shell:  "fastqc -o {params.prefix} --noextract -k 5 -t {threads} -f fastq {input.R1} {input.R2}"

rule cutadapt:
    input: R1 = rawdir + "/"  + "{sample}_R1_001.fastq", R2 = rawdir + "/" + "{sample}_R2_001.fastq"
    output: R1 = "Sample_{sample}/{sample}_trimmed_R1.fastq", R2 = "Sample_{sample}/{sample}_trimmed_R2.fastq"
    params: adapters = adapters
    log: "logs/{sample}_cutadapt.log"
    conda: "envs/cutadapt.yaml"
    threads: 16
    resources: mem_mb=config['medium']['mem'], time=config['medium']['time'], partition=config['medium']['partition']
    shell: "cutadapt -j {threads} -b file:{params.adapters} -B file:{params.adapters} --nextseq-trim=2 --trim-n -n 5 -O 5 -q 10,10 -m 35:35 -o {output.R1} -p {output.R2} {input.R1} {input.R2}"

rule fastqc1:
    input: R1 = "Sample_{sample}/{sample}_trimmed_R1.fastq", R2 = "Sample_{sample}/{sample}_trimmed_R2.fastq"
    output: R1 = "Sample_{sample}/{sample}_trimmed_R1_fastqc.html", R2 = "Sample_{sample}/{sample}_trimmed_R2_fastqc.html"
    params: prefix = "Sample_{sample}/"
    log: "logs/{sample}_fastqc1.log"
    conda: "envs/fastqc.yaml"
    threads: 16
    resources: mem_mb=config['medium']['mem'], time=config['medium']['time'], partition=config['medium']['partition']
    shell: "fastqc -o {params.prefix} --noextract -k 5 -t {threads} -f fastq {input.R1} {input.R2}"

rule seqtk:
    input: R1 = "Sample_{sample}/{sample}_trimmed_R1.fastq", R2 = "Sample_{sample}/{sample}_trimmed_R2.fastq"
    output: R1 = "Sample_{sample}/{sample}_trimmed_R1.sub.fastq", R2 = "Sample_{sample}/{sample}_trimmed_R2.sub.fastq"
    params: n = config['subsamp']
    #log: "logs/{sample}_fastqc1.log"
    conda: "envs/seqtk.yaml"
    #threads: 16
    resources: mem_mb=config['medium']['mem'], time=config['medium']['time'], partition=config['medium']['partition']
    shell: "seqtk sample -s100 {input.R1} {params.n} >{output.R1}; seqtk sample -s100 {input.R2} {params.n} >{output.R2}"

rule fastqscreen:
    input: R1 = "Sample_{sample}/{sample}_trimmed_R1.sub.fastq", R2 = "Sample_{sample}/{sample}_trimmed_R2.sub.fastq"
    output: one = "Sample_{sample}/{sample}_trimmed_R1.sub_screen.png", two = "Sample_{sample}/{sample}_trimmed_R2.sub_screen.png"
    params: prefix = "Sample_{sample}/" , conf = config['fastqscreen']
    conda: "envs/fastqscreen.yaml"
    threads: 16
    resources: mem_mb=config['medium']['mem'], time=config['medium']['time'], partition=config['medium']['partition']
    shell: "fastq_screen --outdir {params.prefix} --threads {threads} --nohits --subset 0  --conf {params.conf} --aligner bowtie2 {input.R1} {input.R2}"

rule kraken2:
    input: R1 = "Sample_{sample}/{sample}_trimmed_R1.sub.fastq", R2 = "Sample_{sample}/{sample}_trimmed_R2.sub.fastq"
    output: report = "Sample_{sample}/{sample}.kraken2.report.txt"
    params: prefix = "Sample_{sample}/{sample}.kraken", kraken2db = config['kraken2db']
    conda: "envs/kraken2.yaml"
    threads: 16
    resources: mem_mb=config['medium']['mem'], time=config['medium']['time'], partition=config['medium']['partition']
    log: "logs/{sample}_kraken.log"
    shell: "kraken2 --fastq-input --threads {threads} --db {params.kraken2db} --output {params.prefix} --report {output.report} --paired {input.R1} {input.R2}"

rule multiqc:
    input: expand("Sample_{sample}/{sample}.sorted.markdup.bam", sample=samples), expand("Sample_{sample}/{sample}.kraken2.report.txt", sample=samples), expand("Sample_{sample}/{sample}_trimmed_R1.sub_screen.png", sample=samples), expand("Sample_{sample}/{sample}_trimmed_R2.sub_screen.png", sample=samples),expand("Genrich/{sample}.frip.score.txt", sample=samples),"Deeptools/bigwig_summary.npz", "ChIPseeker/Distribution_of_Binding_Sites_among_different_ATACseq_data_mqc.png", "Deeptools/plotFingerQualityMetrics.tab", "Deeptools/outRawFragmentLengths.txt"
    output: "multiqc.html"
    conda: "envs/multiqc.yaml"
    resources: mem_mb=config['medium']['mem'], time=config['medium']['time'], partition=config['medium']['partition']
    threads: 16
    shell: "multiqc -f ./ -n {output}"
