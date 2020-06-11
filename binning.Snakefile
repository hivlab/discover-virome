__author__ = "Taavi Päll"
__copyright__ = "Copyright 2020, Avilab"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"

shell.executable("bash")

WRAPPER_PREFIX = "https://raw.githubusercontent.com/avilab/virome-wrappers"

RUNS = ["V554"]

rule all:
    input:
        expand("output/{run}/final.contigs_aln.bam", run = RUNS)

# Fix fasta headers by removing everything after contig id
rule fix_fasta:
    input: 
        "output/{run}/assemble/final.contigs.fa"
    output:
        "output/{run}/contigs-fixed.fa"
    conda:
        "https://raw.githubusercontent.com/avilab/virome-wrappers/master/subset_fasta/environment.yaml"
    script:
        "scripts/fix_fasta.py"

# Map reads to contigs
rule coverage:
    input:
        ref = rules.fix_fasta.output[0], 
        input = "output/{run}/concatenated.fq.gz"
    output:
        out = "output/{run}/final.contigs_aln.bam",
        covstats = "output/{run}/coverage.txt",
        statsfile = "output/{run}/mapcontigs.txt"
    params: 
        extra = lambda wildcards, resources: f"mapper=bbmappacbio maxindel=80 bamscript=bs.sh -Xmx{resources.mem_mb / 1000:.0f}g"
    resources:
        runtime = 1440,
        mem_mb = 8000
    wrapper:
      f"{WRAPPER_PREFIX}/master/bbtools/bbwrap"

