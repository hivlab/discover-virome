
# QC stats
rule fastqc:
    input:
        rules.interleave.output.out
    output:
        html = "output/{run}/fastqc.html",
        zip = "output/{run}/fastqc.zip"
    resources:
        runtime = 120,
        mem_mb = 4000    
    wrapper:
        "0.27.1/bio/fastqc"


rule multiqc:
    input:
        "output/{run}/fastqc.zip",
        "output/{run}/maphost.txt",
        "output/{run}/coverage.txt",
        "output/{run}/mapcontigs.txt"
    output:
        report("output/{run}/multiqc.html", caption = "report/multiqc.rst", category = "Quality control")
    log:
        "output/{run}/log/multiqc.log"
    resources:
        runtime = 120,
        mem_mb = 4000    
    wrapper:
        f"{WRAPPER_PREFIX}/master/multiqc"

