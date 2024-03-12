** Sequencing facility ATAC-Seq pipeline (SF_ATAC-seq) **
The ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) pipeline typically follows several steps to analyze paired-end sequencing data and identify regions of open chromatin (peaks). Here's a basic outline of the process:

## Quality Control (QC):

### [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): Perform initial quality checks on raw sequencing data to assess sequence quality, GC content, over-representation of sequences, etc.
### [Cutadapt] (https://cutadapt.readthedocs.io/en/stable/): Remove low-quality bases, adapter sequences, or other artifacts that may affect downstream analysis.
### Kraken2: Helps detect contamination by identifying unexpected organisms in the sample.
### FastqScreen: Screens sequencing data against a database of known contaminants, such as adapter sequences, PhiX control, and various other sources of contamination. It helps to identify and quantify the presence of these contaminants in the sequencing data.

## Alignment:

### Bowtie2: Aligns paired-end reads to the reference genome.

## Post-alignment Processing:

### Picard MarkDuplicates: Remove duplicate reads introduced during library preparation. 
### Peak Calling: Genrich identifies regions of open chromatin (peaks) 
### Deeptools: assesses the quality of peaks.
### ChIPseeker: identifies the genomic regions associated with open chromatin regions and to perform functional annotation of these regions.

## MultiQC : generates an interactive HTML report that provides a concise summary of the results


