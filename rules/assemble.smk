# --meta   (same as metaspades.py)
#     This flag is recommended when assembling metagenomic data sets
# --only-assembler
#     Runs assembly module only
rule assemble:
    input: 
      rules.cd_hit.output.repres,
      rules.parse_cdhit.output.topn,
      rules.fastq_join.output[0]
    output: 
      "avasta/cdhit/{sample}_cdhit_merged.fa"
      directory("avasta/assemble/{sample}")
    params:
      options = "--meta --only-assembler"
    conda:
      "../envs/spades.yml"
    shell:
      """
      cat {input[0]} {input[1]} > output[0]
	    mkdir -p {output}
	    spades.py {params.options} --merged {input[0]} -s {input[1]} -o {output[1]}
      """
