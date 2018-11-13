
# Split refgenome-filtered reads to joined and unjoined for spades
rule spades_input:
    input:
      rules.fastq_join.output[0],
      rules.refgenomefilter.output.fq
    output:
      join = "avasta/refgenomefilter/{sample}_refgenome_unmapped_join.fq",
      un = "avasta/refgenomefilter/{sample}_refgenome_unmapped_un.fq"
    conda:
      "../envs/biopython.yml"
    script:
      "../scripts/spades_input.py"

# --meta   (same as metaspades.py)
#     This flag is recommended when assembling metagenomic data sets
# --only-assembler
#     Runs assembly module only
rule assemble:
    input: 
      rules.spades_input.output.join,
      rules.spades_input.output.un
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

# Keep only >=500 nt contigs with >=2 coverage
rule filter_contigs:
  input:
    "avasta/assemble/{sample}/scaffolds.fasta"
  output:
    "avasta/fuck/{sample}_good_contigs.fasta"
  params:
    min_length = 500,
    min_coverage = 2
  conda:
    "../envs/biopython.yml"
  script:
    "../scripts/filter_contigs.py"
