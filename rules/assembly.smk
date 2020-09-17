
rule concatenate:
    input:
        rules.merge.output.out, rules.qtrim.output.out
    output:
        temp("output/{group}/{run}/concatenated.fq.gz")
    resources:
        runtime = 120,
        mem_mb = 4000
    shell:
        "cat {input} > {output}"


rule assembly:
    input: 
        se = lambda wildcards: expand("output/{{group}}/{run}/concatenated.fq.gz", run = groups[wildcards.group])
    output: 
        contigs = "output/{group}/assembly/final.contigs.fa"
    params:
        extra = lambda wildcards, resources: f"--presets meta-large --min-contig-len 1000 --verbose -m {resources.mem_mb * 1048576}"
    threads: 8
    log: 
        "output/{group}/log/assembly.log"
    shadow: 
        "minimal"
    resources:
        runtime = lambda wildcards, input: round(609 + 0.03 * input.input.size_mb),
        mem_mb = lambda wildcards, input: round(6000 + 1.5 * input.input.size_mb)
    wrapper:
      f"{WRAPPER_PREFIX}/master/assembly/megahit"


rule fix_fasta:
    input: 
        rules.assembly.output.contigs
    output:
        "output/{group}/contigs-fixed.fa"
    params:
        lambda wildcards: wildcards.group
    conda:
        "https://raw.githubusercontent.com/avilab/virome-wrappers/master/subset_fasta/environment.yaml"
    script:
        "../scripts/fix_fasta.py"


# Calculate assembly coverage stats
# nodisk keeps index in memory, otherwise index will be written once to project root (ref/1) from first run to be processed 
# and reused for other unrelated runs.
# Key "input" will be parsed to "in", "input1" to "in1" etc.
rule mapcontigs:
    input:
        ref = rules.fix_fasta.output[0], 
        input = "output/{group}/{run}/concatenated.fq.gz"
    output:
        out = "output/{group}/{run}/mapcontigs.sam",
        statsfile = "output/{group}/{run}/mapcontigs.txt"
    params: 
        extra = lambda wildcards, resources: f"maxindel=200 strictmaxindel minid=0.9 maxlen=600 nodisk -Xmx{resources.mem_mb / 1000:.0f}g RGLB=lib1 RGPL={PLATFORM} RGID={wildcards.run} RGSM={wildcards.group}"
    resources:
        runtime = 120,
        mem_mb = 40000
    wrapper:
      f"{WRAPPER_PREFIX}/master/bbtools/bbwrap"

rule samtools_sort:
    input:
        rules.mapcontigs.output.out
    output:
        "output/{group}/{run}/sorted.bam"
    params:
        ""
    resources:
        runtime = 120,
        mem_mb = 16000
    threads: 4 # Samtools takes additional threads through its option -@
    wrapper:
        "0.50.4/bio/samtools/sort"


rule samtools_merge:
    input:
        lambda wildcards: expand("output/{{group}}/{run}/sorted.bam", run = groups[wildcards.group])
    output:
        "output/{group}/merged.bam"
    params:
        ""
    resources:
        runtime = 120,
        mem_mb = 16000
    threads:  8  
    wrapper:
        "0.62.0/bio/samtools/merge"


rule samtools_index:
    input:
        "output/{group}/merged.bam"
    output:
        "output/{group}/merged.bam.bai"
    params:
        "" # optional params string
    resources:
        runtime = 120,
        mem_mb = 16000
    wrapper:
        "0.65.0/bio/samtools/index"


rule genomecov:
    input:
        ibam = rules.samtools_merge.output[0]
    output:
        "output/{group}/genomecov.bg"
    params:
        extra = "-bg"
    resources:
        runtime = 120,
        mem_mb = lambda wildcards, attempt: attempt * 8000
    wrapper: 
        f"{WRAPPER_PREFIX}/master/bedtools/genomecov"


# Variant calling
# "Removes any sites with estimated probability of not being polymorphic 
# less than phred 20 (aka 0.01), or probability of polymorphism > 0.99"
# from FreeBayes user manual.
rule lofreq:
    input:
        ref = rules.fix_fasta.output[0],
        bam = rules.samtools_merge.output[0]
    output:
        "output/{group}/lofreq.vcf" 
    params:
        extra="--call-indels --min-cov 50 --max-depth 1000000 --min-bq 30 --min-alt-bq 30 --def-alt-bq 0 --min-mq 20 --max-mq 255 --min-jq 0 --min-alt-jq 0 --def-alt-jq 0 --sig 0.01 --bonf dynamic --no-default-filter"
    resources:
        runtime = 120,
        mem_mb = 4000
    threads: 1
    wrapper:
        f"{WRAPPER_PREFIX}/master/lofreq/call"


rule vcffilter:
    input:
        "output/{group}/lofreq.vcf"
    output:
        "output/{group}/filtered.vcf"
    params:
        extra = "-f 'QUAL > 30 & AF > 0.5'"
    resources:
        runtime = 120,
        mem_mb = 4000
    wrapper:
        f"{WRAPPER_PREFIX}/master/vcflib/vcffilter"
