__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Avilab"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"

# Load libraries
import os
import pandas as pd
from snakemake.utils import validate, makedirs


# Load configuration file with sample and path info
configfile: "config.yaml"
# validate(config, "schemas/config.schema.yaml")


# Load runs and groups
df = pd.read_csv("samples.tsv", sep="\s+", dtype=str).set_index(["sample","run"], drop=False)
validate(df, "schemas/samples.schema.yaml")
samples = df.groupby(level=0).apply(lambda df: df.xs(df.name)["run"].tolist()).to_dict()
SAMPLE = [sample for sample,run in df.index.tolist()]
RUN = [run for sample,run in df.index.tolist()]
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
        expand(["output/{sample}/merged.bam", "output/{sample}/contigs-fixed.fa", "output/{sample}/genomecov.bg", "output/{sample}/filtered.vcf"], sample = list(samples.keys())),
        expand(["output/{sample}/{run}/qtrimmed.fq.gz", "output/{sample}/{run}/fastqc.html"], zip, sample = SAMPLE, run = RUN)
        

def get_fastq(wildcards):
    fq_cols = [col for col in df.columns if "fq" in col]
    fqs = df.loc[(wildcards.sample, wildcards.run), fq_cols].dropna()
    assert len(fq_cols) in [1, 2], "Enter one or two FASTQ file paths"
    if len(fq_cols) == 2:
        return {"in1": fqs[0], "in2": fqs[1]}
    else:
        return {"input": fqs[0]}


include: "rules/preprocess.smk"
include: "rules/assembly.smk"
include: "rules/qc.smk"
