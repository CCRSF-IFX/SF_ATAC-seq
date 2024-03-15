localrules:remove_chrUn, modPeak

rule genrich:
    input: "Sample_{sample}/{sample}.sortedByRead.bam"
    output: genrich = "Genrich/{sample}.narrowPeak.tmp", bedgraph = "Genrich/{sample}.bedgraph", bed = "Genrich/{sample}.bed" 
    params: blacklist = config['blacklist'], prefix = "Genrich/{sample}"
    conda: "envs/genrich.yaml"
    threads: 16
    log: "logs/{sample}.genrich.log"
    resources: mem_mb=config['medium']['mem'], time=config['medium']['time'], partition=config['medium']['partition']
    shell: """
        Genrich -t {input} -o {output.genrich} -j -y -r -v -d 150 -m 5 -e chrM,chrY -E {params.blacklist} -f {params.prefix}.bdg -b {params.prefix}.bed
        cut -f 1,2,3,4 {params.prefix}.bdg | grep -v 'chrUn' | grep -v 'NW' | tail -n +2 > {params.prefix}.bedgraph
        """

rule remove_chrUn:
    input: "Genrich/{sample}.narrowPeak.tmp"
    output: "Genrich/{sample}.narrowPeak"
    shell:"""
        set +e
        grep -v 'chrUn' {input} | grep -v 'NW' > {output}
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            touch {output}
        else
            exit 0
        fi
        """

def file_len(fname):
    with open(fname) as f:
        for i,l in enumerate(f,1):
            pass
    return i


rule NormBdg:
    input: bedgraph = "Genrich/{sample}.bedgraph", bed = "Genrich/{sample}.bed"
    output: normBdg = "Genrich/{sample}.genrich.RPM.bedgraph"
    params: rname = "pl:NormBdg", prefix = "Genrich/{sample}"
    resources: mem_mb=config['medium']['mem'], time=config['medium']['time'], partition=config['medium']['partition']
    threads: 8
    run:
        normFactor = file_len(input.bed)/1e6
        f2 = open(output.normBdg,"w")
        with open(input.bedgraph) as f1:
            for cnt, line in enumerate(f1):
                if 'NA' not in line:
                    line2 = line.strip().split('\t')
                    line2[3] = str( round( float(line2[3]) / normFactor, 6 ) )
                    f2.write( "\t".join(line2) + "\n" )
        f2.close()


rule Bdg2Bw:
    input: bedgraph = "Genrich/{sample}.genrich.RPM.bedgraph"
    output: bigwig = "Genrich/{sample}.genrich.RPM.bw"
    params: sizes = config["sizes"],  prefix = "Genrich/{sample}"
    conda: "envs/genrich.yaml"
    threads: 16
    resources: mem_mb=config['medium']['mem'], time=config['medium']['time'], partition=config['medium']['partition']
    log: "logs/{sample}.Bdg2Bw.log"
    shell: """
            sort -k1,1 -k2,2n -T ./ --parallel={threads} {input.bedgraph} > {params.prefix}.tmp.bedgraph
            bedGraphToBigWig {params.prefix}.tmp.bedgraph {params.sizes} {output.bigwig}
          """


def plotFingerprintInput():
    inputfiles  =  ' '.join(["Sample_%s/%s.sorted.markdup.bam" % (i, i) for i in samples])
    bigwig =  ' '.join(["Genrich/%s.genrich.RPM.bw" % i for i in samples])
    labels = ' '.join(samples)
    return inputfiles, labels, bigwig

rule plotFingerprint:
    input: lambda wildcards: expand("Sample_{sample}/{sample}.sorted.markdup.bam", sample=samples), lambda wildcards: expand("Genrich/{sample}.genrich.RPM.bw", sample=samples)
    output: png = "Deeptools/plotFingerQualityMetrics.png", tab = "Deeptools/plotFingerQualityMetrics.tab"
    params: bam = plotFingerprintInput()[0], labels = plotFingerprintInput()[1], prefix = "Deeptools/plotFingerQualityMetrics"
    conda: "envs/deeptools.yaml"
    threads: 16
    resources: mem_mb=config['medium']['mem'], time=config['medium']['time'], partition=config['medium']['partition']
    log: "logs/plotFingerprint.log"
    shell: """plotFingerprint --numberOfProcessors {threads} -b {params.bam} --labels {params.labels} --minMappingQuality 5 --skipZeros --ignoreDuplicates --outQualityMetrics plotFingerprint_QC_metrics.txt -T "Fingerprints of ATAC-seq samples" --plotFile {output.png} --outQualityMetrics {params.prefix}.tab --outRawCounts Deeptools/RawCounts.tab"""

rule bigwig_summary:
    input: lambda wildcards: expand("Sample_{sample}/{sample}.sorted.markdup.bam", sample=samples), lambda wildcards: expand("Genrich/{sample}.genrich.RPM.bw", sample=samples)
    output: bw = "Deeptools/bigwig_summary.npz"
    params: bigwig = plotFingerprintInput()[2], labels = plotFingerprintInput()[1], bed = config['deeptools_bed']
    conda: "envs/deeptools.yaml"
    log: "logs/bigwig_summary.log"
    threads: 8
    resources: mem_mb=config['medium']['mem'], time=config['medium']['time'], partition=config['medium']['partition']
    shell: "multiBigwigSummary BED-file -p 8 -b {params.bigwig} --labels {params.labels} -o {output.bw} --BED {params.bed}"

rule plotCorrelation:
    input: lambda wildcards: expand("Sample_{sample}/{sample}.sorted.markdup.bam",sample=samples), bw = "Deeptools/bigwig_summary.npz"
    output: spearman = "Deeptools/Spearman_Correlation_of_Samples_Heatmap.png", pearson = "Deeptools/Pearson_Correlation_of_Samples_Heatmap.png",  computeMat = "Deeptools/compute_matrix.gz", heatmap = "Deeptools/Heatmap_of_Gene_Regions.png", profile = "Deeptools/Average_Profile_mqc.png"
    params: bigwig = plotFingerprintInput()[2], labels = plotFingerprintInput()[1], bed = config['deeptools_bed']
    conda: "envs/deeptools.yaml"
    log: "logs/plotCorrelation.log"
    threads: 8
    resources: mem_mb=config['medium']['mem'], time=config['medium']['time'], partition=config['medium']['partition']
    shell: """
        plotCorrelation -in {input.bw} --corMethod spearman --skipZeros --whatToPlot heatmap -T "Spearman Correlation of ATAC-seq samples" -o {output.spearman} --outFileCorMatrix Deeptools/SpearmanCor_bigwigScores.tab
        plotCorrelation -in {input.bw} --corMethod pearson --skipZeros --whatToPlot heatmap -T "Pearson Correlation of ATAC-seq samples" -o {output.pearson} --outFileCorMatrix Deeptools/PearsonCor_bigwigScores.tab
        #computeMatrix scale-regions -p 16 -b 5000 -a 5000 -m 8000 -S {params.bigwig} -R {params.bed} --skipZeros -o {output.computeMat} --samplesLabel {params.labels}
        computeMatrix reference-point -p 16 -a 1000 -b 1000 --referencePoint TSS -S {params.bigwig} -R {params.bed} --skipZeros -o {output.computeMat} --samplesLabel {params.labels}
        plotHeatmap -m {output.computeMat} -out {output.heatmap} --dpi 100
        plotProfile -m {output.computeMat} -out {output.profile} --dpi 100
       """
rule alignmentSieve:
    input: sortBam = "Sample_{sample}/{sample}.sorted.markdup.bam"
    output: bam = "Sample_{sample}/{sample}.sorted.markdup.filtered.bam", sort = "Sample_{sample}/{sample}.sorted.markdup.filtered.sort.bam"
    params: smp = "{sample}"
    conda: "envs/bowtie2.yaml"
    log: "logs/{sample}.alignmentSieve.log"
    threads: 16
    resources: mem_mb=config['medium']['mem'], time=config['medium']['time'], partition=config['medium']['partition']
    shell: """
            alignmentSieve -p {threads} -b {input.sortBam} -o {output.bam} -l {params.smp} --ATACshift --ignoreDuplicates --minMappingQuality 5 --minFragmentLength 20
            samtools sort -@{threads} {output.bam} -o {output.sort}
            samtools index {output.sort}
           """

if len(samples)> 2:
    rule plotPCA:
        input: bw = "Deeptools/bigwig_summary.npz"
        output: "Deeptools/PCA.png"
        shell: "plotPCA -in {input} -o {output} --outFileNameData Deeptools/PCA.tab"


def bamPEFragmentSizeInput():
    inputfiles  =  ' '.join(["Sample_%s/%s.sorted.markdup.filtered.sort.bam" % (i, i) for i in samples])
    labels = ' '.join(samples)
    return inputfiles, labels

rule bamPEFragmentSize:
    input :lambda wildcards: expand("Sample_{sample}/{sample}.sorted.markdup.filtered.sort.bam", sample=samples)
    output: "Deeptools/outRawFragmentLengths.txt"
    params: bamFiles = bamPEFragmentSizeInput()[0], labels = bamPEFragmentSizeInput()[1]
    conda: "envs/deeptools.yaml"
    log: "logs/bamPEFragmentSize.log"
    resources: mem_mb=config['medium']['mem'], time=config['medium']['time'], partition=config['medium']['partition']
    threads: 8
    shell: """bamPEFragmentSize -T "Fragment size of ATAC-seq data" --maxFragmentLength 1000 -b {params.bamFiles} --samplesLabel {params.labels} --outRawFragmentLengths Deeptools/outRawFragmentLengths.txt --table Deeptools/table.txt"""

rule modPeak:
    input: "Genrich/{sample}.narrowPeak"
    output: "Genrich/{sample}.narrowPeak.saf"
    shell: """awk 'BEGIN {{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}} {{print $4, $1, $2+1, $3, "."}}' {input} > {output}"""

rule featureCounts:
    input: saf = "Genrich/{sample}.narrowPeak.saf", bam = "Sample_{sample}/{sample}.sortedByRead.bam"
    output: fc = "Genrich/{sample}.genrich.fc.txt",  err = "Genrich/{sample}.genrich.fc.err",
    conda: "envs/featureCounts.yaml"
    log: "logs/{sample}.featureCounts.log"
    threads: 8
    resources: mem_mb=config['medium']['mem'], time=config['medium']['time'], partition=config['medium']['partition']
    shell:"featureCounts -p -T {threads} -a {input.saf} -F SAF -o {output.fc} {input.bam} 2> {output.err}"

rule fripScore:
    input: "Genrich/{sample}.genrich.fc.err"
    output: "Genrich/{sample}.frip.score.txt"
    run:
        myfile = open(str(output), 'w')
        with open(str(input)) as f:
            for line in f:
                if 'Process BAM file' in line:
                    sample = line.split(' ')[4][:-3]
                if 'Successfully assigned alignments' in line:
                    score = line.split("(")[1].split(')')[0]
                    myfile.write(str(sample) + '\t' + str(score) + '\n')
            myfile.close()

rule ChIPseeker:
    input: expand("Genrich/{sample}.narrowPeak", sample=samples)
    output: "ChIPseeker/Distribution_of_Binding_Sites_among_different_ATACseq_data_mqc.png"
    params: ref = config['reference']
    conda: "envs/R.yaml"
    log: "logs/ChIPseeker.log"
    resources: mem_mb=config['medium']['mem'], time=config['medium']['time'], partition=config['medium']['partition']
    script: "scripts/ChIPseeker.R"

