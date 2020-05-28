
__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Avilab"
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

# Load runs and groups
RUNS = pd.read_csv(config["samples"], sep = "\s+").set_index("run", drop = False)
validate(RUNS, "schemas/samples.schema.yaml")
RUN_IDS = RUNS.index.tolist()
N_FILES = config["split_fasta"]["n_files"]
N = list(range(1, N_FILES + 1, 1))


wildcard_constraints:
    run = "[a-zA-Z0-9]+",
    n = "\d+"

# Main output files
RESULTS = ["viruses.csv", "non-viral.csv", "unassigned.fa"]
BLAST = ["megablast-virus", "blastn-virus", "blastx-virus"] if config["run_blastx"] else ["megablast-virus", "blastn-virus"]
STATS = expand(["output/{run}/multiqc.html"], run = RUN_IDS)
OUTPUTS = expand("output/{run}/{result}", run = RUN_IDS, result = RESULTS) + STATS

# Remote outputs
if config["zenodo"]["deposition_id"]:
    # Load zenodo remote provider module
    from snakemake.remote.zenodo import RemoteProvider as ZENRemoteProvider
    # Setup Zenodo RemoteProvider
    ZEN = ZENRemoteProvider(deposition = config["zenodo"]["deposition_id"], access_token = os.environ["ZENODO_PAT"])
    # Append uploads
    ZENOUTPUTS = ZEN.remote(expand("output/{run}/counts.tgz", run = RUN_IDS))
    OUTPUTS = OUTPUTS + ZENOUTPUTS

# Report
report: "report/workflow.rst"

rule all:
    input:
        OUTPUTS

# Check file exists
def file_exists(file):
    try:
        with open(file, 'r') as fh:
            print("{} is set up correctly".format(file))
    except FileNotFoundError:
        ("Could not find {}").format(file)

# Path to reference genomes
HOST_GENOME = os.getenv("REF_GENOME_HUMAN_MASKED")
# file_exists(REF_GENOME)
TAXON_DB = os.getenv("TAXON_DB")

# Wrappers
WRAPPER_PREFIX = "https://raw.githubusercontent.com/avilab/virome-wrappers/"
BLAST_QUERY = WRAPPER_PREFIX + "master/blast/query"
PARSE_BLAST = WRAPPER_PREFIX + "master/blast/parse"
BLAST_TAXONOMY = WRAPPER_PREFIX + "master/blast/taxonomy"
SUBSET_FASTA = WRAPPER_PREFIX + "master/subset_fasta"

# Rules
include: "rules/preprocess.smk"
include: "rules/blast.smk"
