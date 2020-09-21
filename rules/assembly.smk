
rule concatenate:
    input:
        rules.merge.output.out, rules.qtrim.output.out
    output:
        "output/{group}/{run}/concatenated.fq.gz"
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
        runtime = lambda wildcards, input: round(609 + 0.03 * input.size_mb),
        mem_mb = lambda wildcards, input: round(16000 + 0.07 * input.size_mb ^ 2)
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

