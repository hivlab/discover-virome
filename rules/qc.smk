# QC stats
rule fastqc:
    input:
        rules.interleave.output.out,
    output:
        html="output/{group}/{run}/fastqc.html",
        zip=temp("output/{group}/{run}/fastqc.zip"),
    resources:
        runtime=120,
        mem_mb=4000,
    wrapper:
        "0.27.1/bio/fastqc"


rule multiqc:
    input:
        lambda wildcards: expand(
            [
                "output/{group}/{run}/fastqc.zip",
                "output/{group}/{run}/bhist.txt",
                "output/{group}/{run}/qhist.txt",
                "output/{group}/{run}/aqhist.txt",
                "output/{group}/{run}/bqhist.txt",
                "output/{group}/{run}/lhist.txt",
                "output/{group}/{run}/gchist.txt",
                "output/{group}/{run}/maphost.txt",
            ],
            group=wildcards.group,
            run=[run for group,sample,run in df.index.tolist() if group == wildcards.group],
        ),
    output:
        report(
            "output/{group}/multiqc.html",
            caption="report/multiqc.rst",
            category="Quality control",
        ),
    params:
        "-d -dd 1",
    log:
        "output/{group}/multiqc.log",
    resources:
        runtime=120,
        mem_mb=lambda wildcards, attempt: attempt * 16000,
    wrapper:
        "0.65.0/bio/multiqc"
