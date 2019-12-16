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
from snakemake.utils import validate, makedirs

shell.executable("bash")

# Load configuration file with sample and path info
configfile: "config.yaml"
validate(config, "schemas/config.schema.yaml")

# Load runs and groups
RUNS = pd.read_csv(config["samples"], sep="\s+").set_index("run", drop=False)
validate(RUNS, "schemas/samples.schema.yaml")
RUN_IDS = RUNS.index.tolist()
N_FILES = config["split_fasta"]["n_files"]
N = list(range(1, N_FILES + 1, 1))

# Create slurm logs dir
makedirs("logs/slurm")

wildcard_constraints:
    run = "[a-zA-Z0-9]+",
    n = "\d+",
    blastresult = "[-a-z]+"

# Main output files and target rules
RESULTS = ["viruses.csv", "non-viral.csv", "unassigned.fa"]
BLASTV = ["blastn-virus", "blastx-virus"] if config["run_blastx"] else ["blastn-virus"]
BLASTNR = (
    ["megablast-nt", "blastn-nt", "blastx-nr"]
    if config["run_blastx"]
    else ["megablast-nt", "blastn-nt"]
)
BLAST = BLASTV + BLASTNR
STATS = expand(
    [
        "output/stats/{run}_host-bam-stats.txt",
        "output/stats/{run}_preprocess-stats.tsv",
        "output/stats/{run}_blast.tsv",
        "output/stats/{run}_coverage.txt",
        "output/stats/{run}_basecov.txt",
    ],
    run=RUN_IDS,
)
OUTPUTS = (
    expand(
        ["output/results/{run}_{result}", "output/contigs/{run}_final-contigs.fa"],
        run=RUN_IDS,
        result=RESULTS,
    )
    + STATS
)

# Remote outputs
if config["zenodo"]["deposition_id"]:
    from snakemake.remote.zenodo import RemoteProvider as ZENRemoteProvider
    # Setup Zenodo RemoteProvider
    ZEN = ZENRemoteProvider(deposition = config["zenodo"]["deposition_id"], access_token = os.environ["ZENODO_PAT"])
    ZENOUTPUTS = ZEN.remote(expand(["output/results/{run}_results.tgz", "output/stats/{run}_assembly-stats.tgz", "output/stats/{run}_run-stats.tgz"], run = RUN_IDS))
    OUTPUTS = OUTPUTS + ZENOUTPUTS
    localrules: upload_results, upload_assembly, upload_stats

    rule upload_results:
      input: 
        expand("output/results/{{run}}_{result}", result = RESULTS)
      output: 
        ZEN.remote("output/results/{run}_results.tgz")
      shell: 
        "tar czvf {output} {input}"

    rule upload_stats:
      input: 
        rules.refgenome_bam_stats.output,
        rules.preprocess_stats.output,
        rules.blast_stats.output
      output: 
        ZEN.remote("output/stats/{run}_run-stats.tgz")
      shell: 
        "tar czvf {output} {input}"
    
    rule upload_assembly:
      input:
        rules.assemble_cleanup.output.contigs,
        rules.coverage.output.covstats,
        rules.coverage.output.basecov
      output:
        ZEN.remote("output/stats/{run}_assembly-stats.tgz")
      shell: 
        "tar czvf {output} {input}"

localrules: all
rule all:
    input:
        OUTPUTS

# Path to reference genomes
HOST_GENOME = os.getenv("REF_GENOME_HUMAN")
HOST_TAXID = 9606
TAXON_DB = os.getenv("TAXON_DB")

# Wrappers
wrapper_prefix = "https://raw.githubusercontent.com/avilab/virome-wrappers/"
LN_FILTER = wrapper_prefix + "master/filter/masked"
BWA_UNMAPPED = wrapper_prefix + "master/unmapped"
BLAST_QUERY = wrapper_prefix + "master/blast/query"
PARSE_BLAST = wrapper_prefix + "master/blast/parse"
BLAST_TAXONOMY = wrapper_prefix + "master/blast/taxonomy"
SUBSET_FASTA = wrapper_prefix + "master/subset_fasta"
SEQ_STATS = wrapper_prefix + "master/seqkit/stats"

# Path to Repeatmasker script 
RM = wrapper_prefix + "master/repeatmasker/wrapper.py"

# Modules
include: "rules/preprocess.smk"
include: "rules/blast.smk"

onsuccess:
    email = config["email"]
    shell("mail -s 'Forkflow finished successfully' {email} < {log}")
