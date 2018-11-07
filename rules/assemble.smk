
rule assemble:
    input:
      "avasta/cdhit/{sample}_cdhit.fa",
      "avasta/cdhit/{sample}_cdhit_topn.fa"
    output:

    params:
      outdir: "avasta/assemble"
    conda:
      "../envs/spades.yml"
    script:
      "../scripts/run_spades.py"

