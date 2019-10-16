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

wildcard_constraints:
    run = "[a-zA-Z0-9]+",
    n = "\d+"

# Main output files and target rules
BLAST = ["blastn-virus", "blastx-virus", "megablast-nt", "blastn-nt", "blastx-nr"] if config["run_blastx"] else ["blastn-virus", "megablast-nt", "blastn-nt"]
RESULTS = ["phages-viruses.csv", "non-viral.csv", "query-taxid.csv", "unassigned.fa"]
TAXONOMY = expand("taxonomy/{file}.csv", file = ["names", "nodes", "division"])
STATS = expand(["assemble/stats/{run}_assembly-refgenome-stats.txt", "assemble/stats/{run}_assembly-blast.tsv", "assemble/stats/{run}_coverage.txt", "assemble/stats/{run}_basecov.txt"], run = RUN_IDS)
OUTPUTS = expand(["assemble/results/{run}_{result}", "assemble/contigs/{run}_final-contigs.fa"], run = RUN_IDS, result = RESULTS) + TAXONOMY + STATS

# Remote outputs
if config["zenodo"]["deposition_id"]:
    from snakemake.remote.zenodo import RemoteProvider as ZENRemoteProvider
    # Setup Zenodo RemoteProvider
    ZEN = ZENRemoteProvider(deposition = config["zenodo"]["deposition_id"], access_token = os.environ["ZENODO_PAT"])
    ZENOUTPUTS = ZEN.remote(expand(["assemble/results/{run}_assembly-counts.tgz", "assemble/stats/{run}_assembly-stats.tgz", "assemble/stats/{run}_run-stats.tgz"], run = RUN_IDS))
    OUTPUTS = OUTPUTS + ZENOUTPUTS
    localrules: upload_results, upload_assembly, upload_stats

    rule upload_results:
      input: 
        expand("assemble/results/{{run}}_{result}", result = RESULTS)
      output: 
        ZEN.remote("assemble/results/{run}_assembly-counts.tgz")
      shell: 
        "tar czvf {output} {input}"

    rule upload_stats:
      input: 
        rules.refgenome_bam_stats.output,
        rules.preprocess_stats.output,
        rules.blast_stats.output
      output: 
        ZEN.remote("assemble/stats/{run}_run-stats.tgz")
      shell: 
        "tar czvf {output} {input}"
    
    rule upload_assembly:
      input:
        rules.assemble_cleanup.output.contigs,
        rules.coverage.output.covstats,
        rules.coverage.output.basecov
      output:
        ZEN.remote("assemble/stats/{run}_assembly-stats.tgz")
      shell: 
        "tar czvf {output} {input}"

localrules: all
rule all:
    input:
        OUTPUTS

# Modules
include: "rules/preprocess.smk"
include: "rules/blast.smk"

onsuccess:
    shell("mail -s 'forkflow finished successfully' tapa741@gmail.com < {log}")
