# --meta   (same as metaspades.py)
#     This flag is recommended when assembling metagenomic data sets
# --only-assembler
#     Runs assembly module only
rule assemble:
    input:
        rules.cd_hit.output.repres,
        rules.parse_cdhit.output.topn_fa
    output:
        temp("avasta/cdhit/{sample}_cdhit_merged.fa"),
	directory("avasta/assemble/{sample}")
    params:
        options = "--meta --only-assembler --continue",
	outdir = "avasta/assemble/{sample}"
    conda:
      	"../envs/spades.yml"
    shell:
      	"""
      	cat {input} > {output[0]}
	spades.py {params.options} --merged {ouput[0]} -o {params.outdir}
      	"""

