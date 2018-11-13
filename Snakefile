
__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Avilab"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"

include: "rules/common.smk"

## Target rules
rule all:
    input:
      expand([
      "avasta/assemble/{sample}_good_scaffolds.fasta"
      ], sample = sample_ids)

## Modules
include: "rules/munge.smk"
include: "rules/refgenomefilter.smk"
include: "rules/assemble.smk"
