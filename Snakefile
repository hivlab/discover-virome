__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Avilab"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"

# Load libraries
import os
import pandas as pd
from snakemake.utils import validate, makedirs
from snakemake.utils import min_version

# Check for minimal version
min_version("5.12.3")


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


include: "rules/common.smk"
include: "rules/preprocess.smk"
include: "rules/assembly.smk"
include: "rules/qc.smk"
