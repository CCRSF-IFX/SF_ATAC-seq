# Sequencing Facility ATAC-Seq pipeline (SF_ATAC-seq)

The ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing) pipeline typically follows several steps to analyze paired-end sequencing data and identify regions of open chromatin (peaks). Here's a basic outline of the process:
![SF_ATAC-seq](https://github.com/CCRSF-IFX/SF_ATAC-seq/blob/main/resource/ATAC-Seq.png)

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

## Usage

### Step 1: Obtain a Copy of the Workflow

Clone the Repository: Clone the new repository to your local machine, choosing the directory where you want to perform data analysis. Instructions for cloning can be found [here](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository).

### Step 2: Configure the Workflow

Tailor the workflow to your project's requirements:

Edit `config.yaml` in the `config/` directory to set up the workflow execution parameters.

### Step 3: Load the snakemake version 8 or above

`module load snakemake/8.4.8`

### Step 4: Create a conda environmet

`conda create -n $NAME`

### Step 5: Execute the Workflow

Activate the Conda Environment:

`conda activate $NAME`

Install `mamba`

`conda install -c conda-forge mamba`

Test the Configuration: Perform a dry-run to validate your setup:

`snakemake --use-conda -np`

Local Execution: Execute the workflow on your local machine using $N cores:

`snakemake --use-conda --cores $N`

Here, $N represents the number of cores you wish to allocate for the workflow.

### Cluster Execution: For cluster environments, submit the workflow as follows:

`snakemake --use-conda --cluster slurm --jobs 100`

Replace 100 with the number of jobs you intend to submit simultaneously. Ensure your cluster environment is correctly configured to handle Snakemake jobs.


