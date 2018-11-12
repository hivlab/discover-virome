# --meta   (same as metaspades.py)
#     This flag is recommended when assembling metagenomic data sets
# --only-assembler
#     Runs assembly module only
rule assemble:
    input: 
      rules.cd_hit.output.join,
      rules.cd_hit.output.un
    output: 
      directory("avasta/assemble/{sample}")
    params:
      options = "--meta --only-assembler"
    conda:
      "../envs/spades.yml"
    shell:
      """
	    mkdir -p {output}
	    spades.py {params.options} --merged {input[0]} -s {input[1]} -o {output}
      """
