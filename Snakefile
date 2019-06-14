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
from snakemake.remote.zenodo import RemoteProvider as ZENRemoteProvider
from snakemake.utils import validate

# Load configuration file with sample and path info
configfile: "config.yaml"
validate(config, "schemas/config.schema.yaml")

# Load runs and groups
RUNS = pd.read_csv(config["samples"], sep = "\s+").set_index("run", drop = False)
validate(RUNS, "schemas/samples.schema.yaml")
RUN_IDS = RUNS.index.to_list()
N_FILES = config["split_fasta"]["n_files"]
N = list(range(1, N_FILES + 1, 1))

# Create slurm logs dir
if not os.path.exists("logs/slurm"):
    os.makedirs("logs/slurm")

# Setup Zenodo RemoteProvider
ZEN = ZENRemoteProvider()

wildcard_constraints:
    run = "[a-zA-Z0-9]+",
    n = "\d+"

# Main output files and target rules
RESULTS = ["phages_viruses", "non_viral"]
TAXONOMY = expand("taxonomy/{file}.csv", file = ["names", "nodes", "division"])
STATS = expand(["assemble/stats/{run}_refgenome_stats.txt", "assemble/stats/{run}_blast.tsv", "assemble/stats/{run}_coverage.txt"], run = RUN_IDS)
OUTPUTS = expand(["assemble/results/{run}_query_taxid.csv", "assemble/results/{run}_unassigned.fa", "assemble/results/{run}_{result}.csv", "assemble/{run}/viral.contigs.fa"], run = RUN_IDS, result = RESULTS) + TAXONOMY + STATS

# Remote outputs
if config["zenodo"]["deposition_id"]:
    ZENOUTPUTS = [ZEN.remote(expand("{deposition_id}/files/results/{run}_{result}.csv.tar.gz", deposition_id = config["zenodo"]["deposition_id"], run = RUN_IDS, result = RESULTS)),
    ZEN.remote(expand("{deposition_id}/files/results/{run}_unassigned.fa.tar.gz", deposition_id = config["zenodo"]["deposition_id"], run = RUN_IDS))]

rule all:
    input:
        OUTPUTS, ZENOUTPUTS if config["zenodo"]["deposition_id"] else OUTPUTS

# Modules
include: "rules/preprocess.smk"
include: "rules/blast.smk"
