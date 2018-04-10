
## Use os.path to update file paths from config file
import os.path
import pandas as pd
shell.executable("bash")

## Load configuration file with sample and path info
configfile: "config.yml"
samples = pd.read_table(config["samples"], sep = ",", index_col = "sample", dtype = str)

## Target rule
rule all:
    input:
      expand(os.path.join(config["outdir"], "{sample}/10_repeatmasker_good/{mask}.{n}.fa"), sample = "I1164_12629_Harvard_SIV_196_06_2_24_12_mini", n = [1,2], mask = ["masked", "unmasked"])

include: "rules/munge.smk"
include: "rules/mask.smk"