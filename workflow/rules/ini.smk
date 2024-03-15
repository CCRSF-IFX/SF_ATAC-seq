from snakemake.utils import validate
import pandas as pd
import os
import glob

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
#singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
#validate(config, schema="../schemas/config.schema.yaml")

#validate(samplepd, schema="../schemas/samples.schema.yaml")

analysisdir = config['analysis']
rawdir = config['unaligned']
smpl =  [os.path.basename(file).split('.')[0].split('_R1')[0] for file in glob.glob(rawdir + '/*') if 'R1' in file]
samples = [s.replace('Sample_', '') for s in smpl]
adapters = config["adapters"]
os.chdir(analysisdir)
print (os.getcwd())
print (rawdir)
print (samples)
partition=config['medium']['partition']
print (partition)
