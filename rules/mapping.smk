
# Calculate assembly coverage stats
# nodisk keeps index in memory, otherwise index will be written once to project root (ref/1) from first run to be processed 
# and reused for other unrelated runs.
# Key "input" will be parsed to "in", "input1" to "in1" etc.
rule mapcontigs:
    input:
        ref = rules.fix_fasta.output[0], 
        input = lambda wildcards: expand("output/{{group}}/{run}/concatenated.fq.gz", run = samples[wildcards.sample])
    output:
        out = "output/{group}/{sample}/mapcontigs.bam",
        sorted = "output/{group}/{sample}/mapcontigs_sorted.bam",
        statsfile = "output/{group}/{sample}/mapcontigs.txt"
    params: 
        extra = lambda wildcards, resources: f"maxindel=80 strictmaxindel minid=0.9 maxlen=600 nodisk -Xmx{resources.mem_mb}m bamscript=bs.sh"
    resources:
        runtime = lambda wildcards, attempt: attempt * 240,
        mem_mb = 40000
    wrapper:
      f"{WRAPPER_PREFIX}/master/bbtools/bbwrap"


rule pileup:
    input:
        input = rules.mapcontigs.output.sorted,
        ref = rules.fix_fasta.output[0]
    output:
        out = "output/{group}/{sample}/covstats.txt",
        rpkm = "output/{group}/{sample}/rpkm.txt",
        bincov = "output/{group}/{sample}/bincov.txt"
    params:
        extra = lambda wildcards, resources: f"headerpound=f -Xmx{resources.mem_mb}m"
    resources:
        runtime = lambda wildcards, attempt: attempt * 120,
        mem_mb = lambda wildcards, attempt: attempt * 8000
    wrapper: 
        f"{WRAPPER_PREFIX}/master/bbtools/pileup"


# Variant calling
rule lofreq:
    input:
        ref = rules.fix_fasta.output[0],
        bam = rules.mapcontigs.output.sorted
    output:
        "output/{group}/{sample}/lofreq.vcf" 
    params:
        extra="--call-indels --min-cov 10 --max-depth 1000000 --min-bq 30 --min-alt-bq 30 --def-alt-bq 0 --min-mq 20 --max-mq 255 --min-jq 0 --min-alt-jq 0 --def-alt-jq 0 --sig 0.01 --bonf dynamic --no-default-filter"
    resources:
        runtime = 120,
        mem_mb = 4000
    threads: 1
    wrapper:
        f"{WRAPPER_PREFIX}/master/lofreq/call"
