
# Load required libraries
import os
import json
import glob
import pandas as pd
shell.executable("bash")

# Load configuration file with sample and path info
configfile: "config.yml"
samples = pd.read_table(config["samples"], sep = "\s+", index_col = "sample", dtype = str)
sample_ids = samples.index.values.tolist()
n_files = config["split_fasta"]["n_files"]

# Create slurm logs dir
if not os.path.exists("logs/slurm"):
    os.makedirs("logs/slurm")

