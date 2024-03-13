rule bowtie2:
        input: R1 = "Sample_{sample}/{sample}_trimmed_R1.fastq", R2 = "Sample_{sample}/{sample}_trimmed_R2.fastq"
        output:  final = "Sample_{sample}/{sample}.sorted.bam"
        params: ref = config['ref']
        conda: "envs/bowtie2.yaml" 
        threads: 16
        log: "logs/{sample}_bowtie2.log"
        shell: "bowtie2 -p {threads} -k 10 -x {params.ref} -1 {input.R1} -2 {input.R2} | samtools view -bS - | samtools sort -m 2G -@ {threads} - > {output.final}"

rule sortByRead:
    input: "Sample_{sample}/{sample}.sorted.bam"
    output: final = "Sample_{sample}/{sample}.sortedByRead.bam"
    conda: "envs/bowtie2.yaml"
    threads: 16
    log: "logs/{sample}.sortByRead.log"
    shell: "samtools sort {input} -n -m 2G -@ {threads} -o {output.final}"

rule index:
    input: "Sample_{sample}/{sample}.sorted.bam"
    output: "Sample_{sample}/{sample}.sorted.bam.bai"
    conda: "envs/bowtie2.yaml"
    threads: 16
    log: "logs/{sample}.samtools_index.log"
    shell: "samtools index -@ {threads} {input}"

rule markdup:
    input: "Sample_{sample}/{sample}.sorted.bam"
    output: out = "Sample_{sample}/{sample}.sorted.markdup.bam", metric = "Sample_{sample}/{sample}_MARKEDUPmetrics.txt"
    threads: 16
    log: "logs/{sample}.Markdup.log"
    conda: "envs/markdup.yaml"
    resources: mem_mb=1024
    shell: "picard MarkDuplicates INPUT={input} OUTPUT={output.out} METRICS_FILE={output.metric} ASSUME_SORTED=true CREATE_INDEX=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 VALIDATION_STRINGENCY=LENIENT TMP_DIR=./"

rule markdupIndex:
    input: "Sample_{sample}/{sample}.sorted.markdup.bam"
    output: "Sample_{sample}/{sample}.sorted.markdup.bam.bai"
    conda: "envs/bowtie2.yaml"
    shell: "samtools index {input}"

