#Sequencing facility ATAC-Seq pipeline (SF_ATAC-seq)

The ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) pipeline typically follows several steps to analyze paired-end sequencing data and identify regions of open chromatin (peaks). Here's a basic outline of the process:

## Quality Control (QC):

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): Perform initial quality checks on raw sequencing data to assess sequence quality, GC content, over-representation of sequences, etc.

[Cutadapt](https://cutadapt.readthedocs.io/en/stable/): Remove low-quality bases, adapter sequences, or other artifacts that may affect downstream analysis.

[Kraken2](https://ccb.jhu.edu/software/kraken2/): Helps detect contamination by identifying unexpected organisms in the sample.

[Fastq Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/): Screens sequencing data against a database of known contaminants, such as adapter sequences, PhiX control, and various other sources of contamination. It helps to identify and quantify the presence of these contaminants in the sequencing data.

## Alignment:

[Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml): Aligns paired-end reads to the reference genome.

## Post-alignment Processing:

[Picard MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard): Remove duplicate reads introduced during library preparation. 

Peak Calling: [Genrich](https://github.com/jsh58/Genrich) identifies regions of open chromatin (peaks) 

[Deeptools](https://deeptools.readthedocs.io/en/develop/): assesses the quality of peaks.

[ChIPseeker](https://bioconductor.org/packages/release/bioc/html/ChIPseeker.html): identifies the genomic regions associated with open chromatin regions and to perform functional annotation of these regions.

[MultiQC](https://multiqc.info/) : generates an interactive HTML report that provides a concise summary of the results


