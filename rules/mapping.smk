# Calculate assembly coverage stats
# nodisk keeps index in memory, otherwise index will be written once to project root (ref/1) from first run to be processed
# and reused for other unrelated runs.
# Key "input" will be parsed to "in", "input1" to "in1" etc.
rule mapcontigs:
    input:
        ref=rules.fix_fasta.output[0],
        input=lambda wildcards: expand(
            "output/{{group}}/{run}/concatenated.fq.gz", run=groups[wildcards.group]
        ),
    output:
        out=temp("output/{group}/mapcontigs.bam"),
        statsfile="output/{group}/mapcontigs.txt",
    params:
        extra=(
            lambda wildcards: f"maxlen=600 nodisk RGLB=lib1 RGPL=illumina RGID={wildcards.group} RGSM={wildcards.group}"
        ),
    shadow:
        "minimal"
    resources:
        runtime=lambda wildcards, attempt: attempt * 480,
        mem_mb=40000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.5/bbtools/bbwrap"


rule sort_and_index:
    input:
        rules.mapcontigs.output.out,
    output:
        sorted=temp("output/{group}/mapcontigs_sorted.bam"),
        index=temp("output/{group}/mapcontigs_sorted.bam.bai"),
    params:
        lambda wildcards, resources: f"-m {resources.mem_mb}M",
    threads: 4
    resources:
        mem_mb=16000,
        runtime=lambda wildcards, attempt: attempt * 240,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/samtools/sort_and_index"


rule samtools_index:
    input:
        rules.fix_fasta.output[0],
    output:
        "output/{group}/contigs-fixed.fa.fai",
    params:
        "",
    resources:
        runtime=120,
        mem_mb=4000,
    wrapper:
        "0.68.0/bio/samtools/faidx"


rule lofreq1:
    """
    Variant calling.
    """
    input:
        rules.samtools_index.output[0],
        ref=rules.fix_fasta.output[0],
        bam=rules.sort_and_index.output.sorted,
    output:
        "output/{group}/lofreq1.vcf",
    params:
        extra="",
    resources:
        runtime=120,
        mem_mb=4000,
    threads: 1
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/lofreq/call"


rule indexfeaturefile:
    """
    Index vcf vile.
    """
    input:
        "output/{group}/lofreq1.vcf",
    output:
        "output/{group}/lofreq1.vcf.idx",
    params:
        extra="",
    resources:
        runtime=120,
        mem_mb=4000,
    threads: 1
    wrapper:
        f"{WRAPPER_PREFIX}/v0.5/gatk/indexfeaturefile"


rule sequencedict:
    input:
        rules.fix_fasta.output[0],
    output:
        "output/{group}/contigs-fixed.dict",
    log:
        "output/{group}/log/sequencedict.log",
    resources:
        runtime=120,
        mem_mb=4000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.5/picard/createsequencedictionary"


rule gatk_baserecalibrator:
    input:
        ref=rules.fix_fasta.output[0],
        dict=rules.sequencedict.output[0],
        bam=rules.sort_and_index.output.sorted,
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
    """
    Inserts indel qualities into BAM.
    """
    input:
        ref=rules.fix_fasta.output[0],
        bam=rules.sort_and_index.output.sorted,
        recal_table=rules.gatk_baserecalibrator.output[0],
    output:
        bam="output/{group}/recalibrated.bam",
    resources:
        runtime=120,
        mem_mb=4000,
    wrapper:
        "0.68.0/bio/gatk/applybqsr"


rule lofreq2:
    """
    Variant calling.
    """
    input:
        ref=rules.fix_fasta.output[0],
        bam=rules.applybqsr.output.bam,
    output:
        "output/{group}/lofreq.vcf",
    params:
        extra="--min-cov 10 --max-depth 1000000  --min-bq 30 --min-alt-bq 30 --min-mq 20 --max-mq 255 --min-jq 0 --min-alt-jq 0 --def-alt-jq 0 --sig 0.01 --bonf dynamic --no-default-filter",
    resources:
        runtime=120,
        mem_mb=4000,
    threads: 1
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/lofreq/call"


rule vcffilter:
    """
    Filter variants based on allele frequency.
    """
    input:
        rules.lofreq2.output[0],
    output:
        "output/{group}/filtered.vcf",
    params:
        extra="-f 'AF>0.5'",
    resources:
        runtime=120,
        mem_mb=4000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/vcflib/vcffilter"


rule polish:
    """
    Generate consensus genome, 
    mask positions with low coverage.
    """
    input:
        ref=rules.fix_fasta.output[0],
        bam=rules.applybqsr.output.bam,
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
    """
    Calculate coverage.
    """
    input:
        input=rules.applybqsr.output.bam,
        ref=rules.polish.output.consensus,
    output:
        out="output/{group}/covstats.txt",
        basecov="output/{group}/basecov.txt",
    params:
        extra="concise",
    resources:
        runtime=lambda wildcards, attempt: attempt * 120,
        mem_mb=lambda wildcards, attempt: attempt * 8000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.5/bbtools/pileup"
