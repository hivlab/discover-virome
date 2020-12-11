# Calculate assembly coverage stats
# nodisk keeps index in memory, otherwise index will be written once to project root (ref/1) from first run to be processed
# and reused for other unrelated runs.
# Key "input" will be parsed to "in", "input1" to "in1" etc.
rule mapcontigs:
    input:
        ref=rules.fix_fasta.output[0],
        input=lambda wildcards: expand(
            "output/{{group}}/{run}/concatenated.fq.gz",
            run=sample_runs[wildcards.sample],
        ),
    output:
        out=temp("output/{group}/{sample}/mapcontigs.bam"),
        statsfile="output/{group}/{sample}/mapcontigs.txt",
    log:
        "output/{group}/{sample}/log/mapcontigs.log",
    params:
        extra=(
            lambda wildcards: f"maxlen=600 nodisk RGLB=lib1 RGPL=illumina RGID={wildcards.sample} RGSM={wildcards.group}"
        ),
    shadow:
        "minimal"
    resources:
        runtime=lambda wildcards, attempt: attempt * 480,
        mem_mb=40000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.5/bbtools/bbwrap"


rule samtools_sort:
    input:
        rules.mapcontigs.output.out,
    output:
        temp("output/{group}/{sample}/sorted.bam"),
    log:
        "output/{group}/{sample}/log/samtools_sort.log",
    params:
        extra=lambda wildcards, resources: f"-m {resources.mem_mb}M",
        tmp_dir="/tmp/",
    threads: 8
    resources:
        mem_mb=16000,
        runtime=lambda wildcards, attempt: attempt * 240,
    wrapper:
        "0.68.0/bio/samtools/sort"


rule samtools_merge:
    input:
        lambda wildcards: expand(
            "output/{{group}}/{sample}/sorted.bam",
            sample=group_samples[wildcards.group],
        ),
    output:
        temp("output/{group}/merged.bam"),
    log:
        "output/{group}/log/samtools_merge.log",
    params:
        "-c",
    threads: 8
    resources:
        mem_mb=4000,
        runtime=120,
    wrapper:
        "0.68.0/bio/samtools/merge"


rule lofreq1:
    input:
        rules.samtools_faidx.output[0],
        ref=rules.fix_fasta.output[0],
        bam=rules.samtools_merge.output[0],
    output:
        "output/{group}/lofreq1.vcf",
    log:
        "output/{group}/log/lofreq1.log",
    params:
        extra="--min-cov 10 --max-depth 1000000  --min-bq 30 --min-alt-bq 30 --min-mq 20 --max-mq 255 --min-jq 0 --min-alt-jq 0 --def-alt-jq 0 --sig 0.01 --bonf dynamic --no-default-filter",
    resources:
        runtime=120,
        mem_mb=4000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.6/lofreq/call"


rule indexfeaturefile:
    input:
        "output/{group}/lofreq1.vcf",
    output:
        "output/{group}/lofreq1.vcf.idx",
    log:
        "output/{group}/log/indexfeaturefile.log",
    params:
        extra="",
    resources:
        runtime=120,
        mem_mb=4000,
    threads: 1
    wrapper:
        f"{WRAPPER_PREFIX}/v0.6.1/gatk/indexfeaturefile"


rule gatk_baserecalibrator:
    input:
        ref=rules.fix_fasta.output[0],
        dict=rules.sequencedict.output[0],
        bam=rules.samtools_merge.output[0],
        known=rules.lofreq1.output[0],
        feature_index=rules.indexfeaturefile.output[0],
    output:
        recal_table="output/{group}/recal_table.grp",
    log:
        "output/{group}/log/baserecalibrator.log",
    resources:
        runtime=120,
        mem_mb=4000,
    wrapper:
        "0.68.0/bio/gatk/baserecalibrator"


rule applybqsr:
    input:
        ref=rules.fix_fasta.output[0],
        bam=rules.samtools_merge.output[0],
        recal_table=rules.gatk_baserecalibrator.output[0],
    output:
        bam="output/{group}/recalibrated.bam",
    log:
        "output/{group}/log/applybqsr.log",
    resources:
        runtime=120,
        mem_mb=4000,
    wrapper:
        "0.68.0/bio/gatk/applybqsr"


rule indelqual:
    input:
        ref=rules.fix_fasta.output[0],
        bam=rules.applybqsr.output.bam,
    output:
        "output/{group}/indelqual.bam",
    log:
        "output/{group}/log/indelqual.log",
    params:
        extra="--verbose",
    resources:
        runtime=120,
        mem_mb=4000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.6/lofreq/indelqual"


rule lofreq2:
    input:
        ref=rules.fix_fasta.output[0],
        bam=rules.indelqual.output[0],
    output:
        "output/{group}/lofreq.vcf",
    log:
        "output/{group}/log/lofreq2.log",
    params:
        extra="--call-indels --min-cov 10 --max-depth 1000000  --min-bq 30 --min-alt-bq 30 --min-mq 20 --max-mq 255 --min-jq 0 --min-alt-jq 0 --def-alt-jq 0 --sig 0.01 --bonf dynamic --no-default-filter",
    resources:
        runtime=120,
        mem_mb=4000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.6/lofreq/call"


rule vcffilter:
    input:
        rules.lofreq2.output[0],
    output:
        "output/{group}/filtered.vcf",
    params:
        extra="-f 'AF > 0.5'",
    resources:
        runtime=120,
        mem_mb=4000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/vcflib/vcffilter"


rule polish:
    input:
        ref=rules.fix_fasta.output[0],
        bam=rules.indelqual.output[0],
        vcf=rules.vcffilter.output[0],
    output:
        vcfgz="output/{group}/filtered.vcf.gz",
        consensus="output/{group}/contigs-polished.fa",
        consensus_masked=temp("output/{group}/contigs-polished-masked.fa"),
        bed=temp("output/{group}/contigs-mask.bed"),
    log:
        "output/{group}/log/polish.log",
    params:
        mask=1,
    resources:
        runtime=120,
        mem_mb=4000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.3/genome-consensus"


rule pileup:
    input:
        input=rules.samtools_sort.output[0],
        ref=rules.fix_fasta.output[0],
    output:
        out="output/{group}/{sample}/covstats.txt",
        basecov="output/{group}/{sample}/basecov.txt",
    log:
        "output/{group}/{sample}/log/pileup.log",
    params:
        extra="concise",
    resources:
        runtime=lambda wildcards, attempt: attempt * 120,
        mem_mb=lambda wildcards, attempt: attempt * 8000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.5/bbtools/pileup"
