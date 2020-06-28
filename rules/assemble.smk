
rule concatenate:
    input:
        rules.merge.output.out, rules.qtrim.output.out
    output:
        temp("output/{run}/concatenated.fq.gz")
    resources:
        runtime = 120,
        mem_mb = 4000
    shell:
        "cat {input} > {output}"


rule assemble:
    input: 
        se = expand("output/{run}/concatenated.fq.gz", run = RUN_IDS)
    output: 
        contigs = "output/assemble/final.contigs.fa"
    params:
        extra = lambda wildcards, resources: f"--min-contig-len 1000 -m {resources.mem_mb * 1048576}"
    threads: 8
    log: 
        "output/assemble/log/assemble.log"
    shadow: 
        "minimal"
    resources:
        runtime = 1440,
        mem_mb = 96000
    wrapper:
      f"{WRAPPER_PREFIX}/master/assembly/megahit"


# Calculate assembly coverage stats
# nodisk keeps index in memory, otherwise index will be written once to project root (ref/1) from first run to be processed 
# and reused for other unrelated runs.
# Key "input" will be parsed to "in", "input1" to "in1" etc.
rule coverage:
    input:
        ref = rules.assemble.output.contigs, 
        input = expand("output/{run}/concatenated.fq.gz", run = RUN_IDS)
    output:
        out = "output/assemble/final.contigs_aln.sam",
        covstats = "output/assemble/coverage.txt",
        statsfile = "output/assemble/mapcontigs.txt"
    params: 
        extra = lambda wildcards, resources: f"mapper=bbmappacbio maxindel=80 strictmaxindel minid=0.9 -Xmx{resources.mem_mb / 1000:.0f}g"
    resources:
        runtime = 1440,
        mem_mb = 16000
    wrapper:
      f"{WRAPPER_PREFIX}/master/bbtools/bbwrap"
