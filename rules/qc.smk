
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
rule bamstats:
    input:
        rules.samtools_merge.output[0]
    output:
        "output/{sample}/bamstats.txt"
    resources:
        runtime = 120,
        mem_mb = 8000
    wrapper:
        "0.42.0/bio/samtools/stats"


rule multiqc:
    input:
        expand(["output/{sample}/{run}/fastqc.zip",], zip, sample = SAMPLE, run = RUN),
        expand(["output/{sample}/bamstats.txt"], sample = samples.keys())
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
      f"{WRAPPER_PREFIX}/master/multiqc"
