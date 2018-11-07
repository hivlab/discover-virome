# --meta   (same as metaspades.py)
#     This flag is recommended when assembling metagenomic data sets
# --only-assembler
#     Runs assembly module only
rule assemble:
    input:
        "avasta/cdhit/{sample}_cdhit.fa",
        "avasta/cdhit/{sample}_cdhit_topn.fa"
    output:
        directory("avasta/assemble/{sample}")
    params:
        options = "--meta --only-assembler --continue"
    conda:
      "../envs/spades.yml"
    shell:
      """
      cat {input} | spades.py {params.options} --merged - -o {output}
      """

