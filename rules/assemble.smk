# --meta   (same as metaspades.py)
#     This flag is recommended when assembling metagenomic data sets
# --only-assembler
#     Runs assembly module only
rule assemble:
    input:
        rules.cd_hit.output.repres,
        rules.parse_cdhit.output.topn_fa
    output:
    	merged = temp("avasta/cdhit/{sample}_cdhit_merged.fa"),
	outdir = directory("avasta/assemble/{sample}")
    params:
        options = "--meta --only-assembler"
    conda:
      	"../envs/spades.yml"
    shell:
      	"""
      	cat {input} > {output.merged}
	mkdir -p {output.outdir}
	spades.py {params.options} --merged {output.merged} -o {output.outdir}
      	"""

