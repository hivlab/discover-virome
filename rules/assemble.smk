# --meta   (same as metaspades.py)
#     This flag is recommended when assembling metagenomic data sets
# --only-assembler
#     Runs assembly module only
rule assemble:
    input: 
      rules.refgenomefilter.output.fq,
      rules.fastq_join.output[1]
    output: directory("avasta/assemble/{sample}")
    params:
      options = "--meta --only-assembler"
    conda:
      "../envs/spades.yml"
    shell:
      """
	    mkdir -p {output}
	    spades.py {params.options} --merged {input[0]} -s {input[1]} -o {output}
      """
