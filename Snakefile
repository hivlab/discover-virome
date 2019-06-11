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

wildcard_constraints:
    run = "[a-zA-Z0-9]+"

# Main output files and target rules
RESULTS = ["phages_viruses", "non_viral"]
TAXONOMY = expand("taxonomy/{file}.csv", file = ["names", "nodes", "division"])
STATS = expand(["assemble/stats/{run}_refgenome_stats.txt", "assemble/stats/{run}_blast.tsv", "assemble/stats/{run}_coverage.txt"], run = RUN_IDS)
OUTPUTS = expand(["assemble/results/{run}_query_taxid.csv", "assemble/results/{run}_unassigned.fa", "assemble/results/{run}_{result}.csv"], run = RUN_IDS, result = RESULTS) + TAXONOMY + STATS

rule all:
    input:
        OUTPUTS, 
        ZEN.remote(expand(["{deposition_id}/files/assemble/results/{run}_query_taxid.csv", "{deposition_id}/files/assemble/results/{{run, [^_]+}}_{{result}}.{{ext}}"], run = RUN_IDS, result = RESULTS, deposition_id = config["zenodo"]["deposition_id"])) if config["zenodo"]["deposition_id"] else OUTPUTS

# Modules
include: "rules/preprocess.smk"
include: "rules/blast.smk"
