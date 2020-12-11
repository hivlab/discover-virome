
rule concatenate:
    input:
        rules.merge.output.out,
        rules.qtrim.output.out,
    output:
        temp("output/{group}/{run}/concatenated.fq.gz"),
    resources:
        runtime=120,
        mem_mb=4000,
    shell:
        "cat {input} > {output}"


rule assembly:
    input:
        se=lambda wildcards: expand(
            "output/{{group}}/{run}/concatenated.fq.gz", run=group_runs[wildcards.group]
        ),
    output:
        contigs=temp("output/{group}/assembly/final.contigs.fa"),
    log:
        "output/{group}/log/assembly.log",
    params:
        extra=(
            lambda wildcards, resources: f"--presets meta-large --min-contig-len 1000 --verbose -m {resources.mem_mb * 1048576}"
        ),
    shadow:
        "minimal"
    resources:
        runtime=lambda wildcards, input: round(2600 + 0.06 * input.size_mb),
        mem_mb=lambda wildcards, input: round(20000 + 2.22 * input.size_mb),
    threads: 8
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/assembly/megahit"


rule fix_fasta:
    input:
        rules.assembly.output.contigs,
    output:
        "output/{group}/contigs-fixed.fa",
    params:
        lambda wildcards: wildcards.group,
    conda:
        f"{WRAPPER_PREFIX}/v0.2/subset_fasta/environment.yaml"
    script:
        "../scripts/fix_fasta.py"


rule samtools_faidx:
    input:
        rules.fix_fasta.output[0],
    output:
        "output/{group}/contigs-fixed.fa.fai",
    log:
        "output/{group}/log/samtools_faidx.log",
    params:
        "",
    resources:
        runtime=120,
        mem_mb=4000,
    wrapper:
        "0.68.0/bio/samtools/faidx"


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
