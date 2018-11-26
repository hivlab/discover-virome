
## Split reads to smaller chunks for Repeatmasker
rule split_fasta:
  input:
    rules.cd_hit.output.repres
  output:
    expand("avasta/mask/{{sample}}_repeatmasker_{n}.fa", n = list(range(1, n_files + 1, 1)))
  params:
    config["split_fasta"]["n_files"]
  conda:
    "../envs/biopython.yaml"
  script:
    "../scripts/split_fasta.py"

os.environ['REPEATMASKER_REPBASE_FILE']=config["repbase_file"]

## Repeatmasker
# Outputs are generated from input file names by RepeatMasker
# must have file extension '.masked'
# If no repetitive sequences were detected symlink output to input file
rule repeatmasker:
  input:
    "avasta/mask/{sample}_repeatmasker_{n}.fa"
  output:
    masked = "avasta/mask/{sample}_repeatmasker_{n}.fa.masked",
    out = "avasta/mask/{sample}_repeatmasker_{n}.fa.out",
    tbl = "avasta/mask/{sample}_repeatmasker_{n}.fa.tbl"
  params:
    outdir = "avasta/mask"
  threads: 8
  conda: "../envs/repeatmasker.yaml"
  shell:
    """
    RepeatMasker -qq -pa {threads} {input} -dir {params.outdir}
    if head -n 1 {output.out} | grep -q "There were no repetitive sequences detected"
      then ln -sr {input} {output.masked}
    fi
    """

## Filter repeatmasker output [10]
# 1) Sequences > 50 nt of consecutive sequence without N
# 2) Sequences with >= 40% of total length of being masked
# input, output, and params names must match function arguments
rule repeatmasker_good:
  input:
    masked = rules.repeatmasker.output.masked,
    original = rules.repeatmasker.input
  output:
    masked_filt = "avasta/mask/{sample}_repmaskedgood_{n}.fa",
    original_filt = "avasta/mask/{sample}_unmaskedgood_{n}.fa"
  params:
    min_length = 50,
    por_n = 40
  conda:
    "../envs/biopython.yaml"
  script:
    "../scripts/filter_masked.py"
