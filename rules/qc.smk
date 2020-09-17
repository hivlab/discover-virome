
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


rule multiqc:
    input:
        expand(["output/{group}/{run}/bhist.txt", "output/{group}/{run}/qhist.txt", "output/{group}/{run}/aqhist.txt", "output/{group}/{run}/bqhist.txt", "output/{group}/{run}/lhist.txt", "output/{group}/{run}/gchist.txt", "output/{group}/{run}/maphost.txt"], zip, group = GROUP, run = RUN)
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
