
__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Avilab"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"

# Load libraries
import os
import json
import glob
import pandas as pd
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.utils import validate
shell.executable("bash")

# Load configuration file with sample and path info
configfile: "config.yaml"
validate(config, "schemas/config.schema.yaml")
SAMPLES = pd.read_table(config["samples"], sep = "\s+").set_index("sample", drop=False)
validate(SAMPLES, "schemas/samples.schema.yaml")
SAMPLE_IDS = SAMPLES.index.values.tolist()
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

# Create slurm logs dir
if not os.path.exists("logs/slurm"):
    os.makedirs("logs/slurm")

rule all:
    input:
        expand(["assemble/{sample}/final.contigs.fa", 
                "align/{sample}/aln.sam.gz", 
                "align/{sample}/coverage.tsv",
                "align/{sample}_sorted.bam",
                "network/{sample}/network.txt"], sample = SAMPLE_IDS)

# Modules
include: "rules/trim.smk"
include: "rules/assemble.smk"
include: "rules/network.smk"
