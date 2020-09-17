__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Avilab"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"

# Load libraries
import os
import pandas as pd
from snakemake.utils import validate, makedirs


# Load configuration file with group and path info
configfile: "config.yaml"
# validate(config, "schemas/config.schema.yaml")


# Load samples
df = pd.read_csv("samples.tsv", sep="\s+", dtype=str).set_index(["group","run"], drop=False)
validate(df, "schemas/samples.schema.yaml")
groups = df.groupby(level=0).apply(lambda df: df.xs(df.name)["run"].tolist()).to_dict()
GROUP = [group for group,run in df.index.tolist()]
RUN = [run for group,run in df.index.tolist()]
PLATFORM = config["platform"]

# Path to reference genomes
HOST_GENOME = os.getenv("REF_GENOME_HUMAN_MASKED")
TAXON_DB = os.getenv("TAXON_DB")


# Wrappers
# Wrappers repo: https://github.com/avilab/virome-wrappers
WRAPPER_PREFIX = "https://raw.githubusercontent.com/avilab/virome-wrappers"


# Report
report: "report/workflow.rst"


onsuccess:
    email = config["password"]
    shell("mail -s 'Forkflow finished successfully' {email} < {log}")


rule all:
    input: 
        "output/multiqc.html",
        expand(["output/{group}/contigs-fixed.fa"], group = list(groups.keys())),
        expand(["output/{group}/{run}/fastqc.html"], zip, group = GROUP, run = RUN)
        

def get_fastq(wildcards):
    fq_cols = [col for col in df.columns if "fq" in col]
    fqs = df.loc[(wildcards.group, wildcards.run), fq_cols].dropna()
    assert len(fq_cols) in [1, 2], "Enter one or two FASTQ file paths"
    if len(fq_cols) == 2:
        return {"in1": fqs[0], "in2": fqs[1]}
    else:
        return {"input": fqs[0]}


include: "rules/preprocess.smk"
include: "rules/assembly.smk"
include: "rules/qc.smk"
