
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
        se = rules.concatenate.output[0]
    output: 
        contigs = "output/{run}/assemble/final.contigs.fa"
    params:
        extra = "--min-contig-len 1000"
    threads: 8
    log: 
        "output/{run}/log/assemble.log"
    shadow: 
        "minimal"
    resources:
        runtime = lambda wildcards, attempt: attempt * 120,
        mem_mb = 8000
    wrapper:
      f"{WRAPPER_PREFIX}/master/assembly/megahit"


# Calculate assembly coverage stats
# nodisk keeps index in memory, otherwise index will be written once to project root (ref/1) from first run to be processed 
# and reused for other unrelated runs.
# Key "input" will be parsed to "in", "input1" to "in1" etc.
rule coverage:
    input:
        ref = rules.assemble.output.contigs, 
        input = rules.concatenate.output[0]
    output:
        out = "output/{run}/final.contigs_aln.sam",
        covstats = "output/{run}/coverage.txt",
        statsfile = "output/{run}/mapcontigs.txt"
    params: 
        extra = lambda wildcards, resources: f"mapper=bbmappacbio maxindel=80 nodisk -Xmx{resources.mem_mb / 1000:.0f}g"
    resources:
        runtime = 1440,
        mem_mb = 4000
    wrapper:
      f"{WRAPPER_PREFIX}/master/bbtools/bbwrap"
