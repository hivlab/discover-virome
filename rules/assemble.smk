
rule assemble:
    input:
        "avasta/cdhit/{sample}_cdhit.fa",
        "avasta/cdhit/{sample}_cdhit_topn.fa"
    output:
        directory("avasta/assemble/{sample}")
    params:
        options = "--meta --only-assembler"
    conda:
      "../envs/spades.yml"
    shell:
      """
      spades.py {params.options} -o {output}
      """

