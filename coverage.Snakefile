
samples = {

    "A": ["AL1", "AL2"]

}

WRAPPER_PREFIX = "https://raw.githubusercontent.com/avilab/virome-wrappers"

rule all:
    input: expand(["output/merged-{sample}/contigs.bg", "output/merged-{sample}/contigs.rpkm"], sample = samples.keys())

rule samtools_merge:
    input:
        lambda wildcards: expand("output/{{sample}}/{run}/contigs_sorted.bam", run = samples[wildcards.sample])
    output:
        "output/merged-{sample}/contigs.bam"
    params:
        ""
    threads:  8 
    resources:
        runtime = 120,
        mem_mb = lambda wildcards, attempt: attempt * 8000
    wrapper:
        "0.62.0/bio/samtools/merge"


rule pileup:
    input:
        input = rules.samtools_merge.output[0],
        ref = "output/assemble/contigs-fixed.fa"
    output:
        out = "output/merged-{sample}/contigs.bg",
        fpkm = "output/merged-{sample}/contigs.rpkm"
    params:
        extra = lambda wildcards, resources: f"-Xmx{resources.mem_mb / 1000:.0f}g"
    resources:
        runtime = 120,
        mem_mb = lambda wildcards, attempt: attempt * 8000
    wrapper: 
        f"{WRAPPER_PREFIX}/master/bbtools/pileup"
