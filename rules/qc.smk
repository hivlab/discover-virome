
# QC stats
rule fastqc:
    input:
        unpack(get_fastq)
    output:
        html = "output/{sample}/{run}/fastqc.html",
        zip = "output/{sample}/{run}/fastqc.zip"
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
        "output/{sample}/samtools-stats.txt"
    resources:
        runtime = 120,
        mem_mb = 8000
    wrapper:
        "0.42.0/bio/samtools/stats"


rule samtools_flagstat:
    input:
        rules.samtools_merge.output[0]
    output:
        "output/{sample}/samtools-flagstats.txt"
    wrapper:
        "0.65.0/bio/samtools/flagstat"


rule samtools_idxstats:
    input:
        bam = rules.samtools_merge.output[0],
        idx = rules.samtools_index.output[0]
    output:
        "output/{sample}/samtools-idxstats.txt"
    log:
        "output/{sample}/log/idxstats.log"
    wrapper:
        "0.65.0/bio/samtools/idxstats"


rule multiqc:
    input:
        expand(["output/{sample}/samtools-stats.txt", "output/{sample}/samtools-flagstats.txt", "output/{sample}/samtools-idxstats.txt"], sample = samples.keys())
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
