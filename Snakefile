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


validate(config, "schemas/config.schema.yaml")

# Load samples
def run_index_dict(df, level=0, slice="run"):
    L = [{k: v.to_list()} for k, v in list(df.groupby(level=level)[slice])]
    return {k: v for d in L for k, v in d.items()}


df = (
    pd.read_csv("samples.tsv", sep="\s+", dtype=str)
    .set_index(["group", "sample", "run"], drop=False)
    .sort_index()
)
validate(df, "schemas/samples.schema.yaml")

# Sequencing run groups for contigs co-assembly
group_runs = run_index_dict(df, level=0)

# Sequencing run groups/biosamples for contigs mapping
sample_runs = run_index_dict(df, level=1)

# Samples in group
group_samples = run_index_dict(df, level=0, slice="sample")

GROUP = [group for group, sample, run in df.index.tolist()]
SAMPLE = [sample for group, sample, run in df.index.tolist()]
RUN = [run for group, sample, run in df.index.tolist()]
N = list(range(0, config["splits"], 1))

# Path to reference genomes
HOST_GENOME = os.getenv("REF_GENOME_HUMAN_MASKED")
TAXON_DB = os.getenv("TAXON_DB")

# Wrappers Github repo: https://github.com/avilab/virome-wrappers
WRAPPER_PREFIX = "https://raw.githubusercontent.com/avilab/virome-wrappers"


# Report
report: "report/workflow.rst"


onsuccess:
    email = config["email"]
    shell("mail -s 'Forkflow finished successfully' {email} < {log}")


wildcard_constraints:
    group="[^/]",


rule all:
    input:
        expand(
            [
                "output/{group}/lofreq.vcf",
                "output/{group}/multiqc.html",
                "output/{group}/contigs-fixed.fa",
                "output/{group}/viruses.csv",
                "output/{group}/non-viral.csv",
                "output/{group}/unassigned.fa",
            ],
            group=list(group_runs.keys()),
        ),
        expand(
            ["output/{group}/{sample}/covstats.txt"], zip, group=GROUP, sample=SAMPLE
        ),


include: "rules/common.smk"
include: "rules/preprocess.smk"
include: "rules/assembly.smk"
include: "rules/mapping.smk"
include: "rules/qc.smk"
include: "rules/mask.smk"
include: "rules/blast.smk"
