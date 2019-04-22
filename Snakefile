__author__ = "Taavi Päll"
__copyright__ = "Copyright 2019, Avilab"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"

# Load libraries
import os
import json
import glob
import pandas as pd
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.utils import validate

# Load configuration file with sample and path info
configfile: "config.yaml"
#validate(config, "schemas/config.schema.yaml")
SAMPLES = pd.read_csv(config["samples"], sep = "\\s+").set_index("run", drop = False)
#validate(SAMPLES, "schemas/samples.schema.yaml")
SAMPLE_IDS, RUN_IDS = list(SAMPLES["sample"]), SAMPLES.run.to_list()
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

# Create slurm logs dir
if not os.path.exists("logs/slurm"):
    os.makedirs("logs/slurm")

rule all:
    input:
        expand(["assemble/{run}/final.contigs.fa", 
        "assemble/{run}/coverage.txt",
        "mask/{run}_repmaskedgood.fa",
        "mask/{run}_unmaskedgood.fa"], run = RUN_IDS)

# Modules
include: "rules/trim.smk"
include: "rules/assemble.smk"
include: "rules/mask.smk"

