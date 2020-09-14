
rule concatenate:
    input:
        rules.merge.output.out, rules.qtrim.output.out
    output:
        temp("output/{sample}/{run}/concatenated.fq.gz")
    resources:
        runtime = 120,
        mem_mb = 4000
    shell:
        "cat {input} > {output}"


rule assembly:
    input: 
        se = lambda wildcards: expand("output/{{sample}}/{run}/concatenated.fq.gz", run = samples[wildcards.sample])
    output: 
        contigs = "output/{sample}/assembly/final.contigs.fa"
    params:
        extra = lambda wildcards, resources: f"--min-contig-len 1000 -m {resources.mem_mb * 1048576}"
    threads: 8
    log: 
        "output/{sample}/log/assembly.log"
    shadow: 
        "minimal"
    resources:
        runtime = lambda wildcards, attempt: attempt * 600,
        mem_mb = 36000
    wrapper:
      f"{WRAPPER_PREFIX}/master/assembly/megahit"


rule fix_fasta:
    input: 
        rules.assembly.output.contigs
    output:
        "output/{sample}/contigs-fixed.fa"
    params:
        lambda wildcards: wildcards.sample
    conda:
        "https://raw.githubusercontent.com/avilab/virome-wrappers/master/subset_fasta/environment.yaml"
    script:
        "../scripts/fix_fasta.py"


checkpoint concatcontigs:
    input:
        lambda wildcards: expand("output/{sample}/contigs-fixed.fa", run = samples[wildcards.sample])
    output:
        temp("output/concatcontigs.fa")
    resources:
        runtime = 120,
        mem_mb = 4000
    shell:
        "cat {input} > {output}"


# Run cd-hit to find and cluster duplicate sequences
rule cd_hit:
    input:
        rules.concatcontigs.output[0]
    output:
        repres = temp("output/cdhit.fa"),
        clstr = temp("output/cdhit.fa.clstr")
    params:
        extra = "-c 0.984 -G 0 -n 10 -d 0 -aS 0.984 -r 1 -M 0"
    log:
        "output/log/cdhit.log"
    threads: 4
    resources:
        runtime = 2880,
        mem_mb = 14000
    wrapper:
        f"{WRAPPER_PREFIX}/master/cdhit"


# Calculate assembly coverage stats
# nodisk keeps index in memory, otherwise index will be written once to project root (ref/1) from first run to be processed 
# and reused for other unrelated runs.
# Key "input" will be parsed to "in", "input1" to "in1" etc.
rule mapcontigs:
    input:
        ref = rules.fix_fasta.output[0], 
        input = "output/{sample}/{run}/concatenated.fq.gz"
    output:
        out = "output/{sample}/{run}/mapcontigs.sam",
        statsfile = "output/{sample}/{run}/mapcontigs.txt"
    params: 
        extra = lambda wildcards, resources: f"maxindel=200 strictmaxindel minid=0.9 maxlen=600 nodisk -Xmx{resources.mem_mb / 1000:.0f}g RGLB=lib1 RGPL={PLATFORM} RGID={wildcards.run} RGSM={wildcards.sample}"
    resources:
        runtime = 120,
        mem_mb = 4000
    wrapper:
      f"{WRAPPER_PREFIX}/master/bbtools/bbwrap"

rule samtools_sort:
    input:
        rules.mapcontigs.output.out
    output:
        "output/{sample}/{run}/sorted.bam"
    params:
        ""
    resources:
        runtime = 120,
        mem_mb = 4000
    threads: 4 # Samtools takes additional threads through its option -@
    wrapper:
        "0.50.4/bio/samtools/sort"


rule samtools_merge:
    input:
        lambda wildcards: expand("output/{{sample}}/{run}/sorted.bam", run = samples[wildcards.sample])
    output:
        "output/{sample}/merged.bam"
    params:
        ""
    threads:  8  
    wrapper:
        "0.62.0/bio/samtools/merge"


rule genomecov:
    input:
        ibam = rules.samtools_merge.output[0]
    output:
        "output/{sample}/genomecov.bg"
    params:
        extra = "-bg"
    resources:
        runtime = 120,
        mem_mb = lambda wildcards, attempt: attempt * 8000
    wrapper: 
        f"{WRAPPER_PREFIX}/master/bedtools/genomecov"
