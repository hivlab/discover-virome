
# QC stats
rule fastqc:
    input:
        unpack(get_fastq)
    output:
        html = "output/{group}/{run}/fastqc.html",
        zip = "output/{group}/{run}/fastqc.zip"
    resources:
        runtime = 120,
        mem_mb = 4000    
    wrapper:
        "0.27.1/bio/fastqc"


# Host mapping stats
rule samtools_stats:
    input:
        rules.samtools_merge.output[0]
    output:
        "output/{group}/samtools-stats.txt"
    resources:
        runtime = 120,
        mem_mb = 8000
    wrapper:
        "0.42.0/bio/samtools/stats"


rule samtools_flagstat:
    input:
        rules.samtools_merge.output[0]
    output:
        "output/{group}/samtools-flagstats.txt"
    wrapper:
        "0.65.0/bio/samtools/flagstat"


rule samtools_idxstats:
    input:
        bam = rules.samtools_merge.output[0],
        idx = rules.samtools_index.output[0]
    output:
        "output/{group}/samtools-idxstats.txt"
    log:
        "output/{group}/log/idxstats.log"
    wrapper:
        "0.65.0/bio/samtools/idxstats"


rule multiqc:
    input:
        expand(["output/{group}/samtools-stats.txt", "output/{group}/samtools-flagstats.txt", "output/{group}/samtools-idxstats.txt"], group = groups.keys())
    output:
        report("output/multiqc.html", caption = "report/multiqc.rst", category = "Quality control")
    params:
        "-d -dd 1"
    log:
        "output/multiqc.log"
    resources:
        runtime = 120,
        mem_mb = 4000    
    wrapper:
        "0.65.0/bio/multiqc"
