
# Tantan mask of low complexity DNA sequences
rule tantan:
  input: rules.assemble.output.contigs
  output:
    "mask/{run}_tantan.fasta"
  conda:
      "../envs/tantan.yaml"
  shell:
    """
    tantan -x N {input} > {output}
    """

# Filter tantan output
# 1) Sequences > 50 nt of consecutive sequence without N
# 2) Sequences with >= 40% of total length of being masked
rule tantan_good:
  input:
    masked = rules.tantan.output
  output:
    masked_filt = "mask/{run}_tantangood.fasta"
  params:
    min_length = 50,
    por_n = 40
  wrapper:
    "https://raw.githubusercontent.com/avilab/snakemake-wrappers/master/filter/masked"


os.environ['REPEATMASKER_REPBASE_FILE']=config["repbase_file"]
# Repeatmasker
# Outputs are generated from input file names by RepeatMasker
# must have file extension '.masked'
# If no repetitive sequences were detected syamlink output to input file
rule repeatmasker:
  input:
    fa = rules.tantan_good.output
  output:
    masked = "mask/{run}_repeatmasker.fa.masked",
    out = "mask/{run}_repeatmasker.fa.out",
    tbl = "mask/{run}_repeatmasker.fa.tbl"
  params:
    outdir = "mask"
  threads: 8
  conda: "../envs/repeatmasker.yaml"
  shell:
    """
    RepeatMasker -qq -pa {threads} {input.fa} -dir {params.outdir}
    if head -n 1 {output.out} | grep -q "There were no repetitive sequences detected"
      then ln -sr {input.fa} {output.masked} \
           && touch {output.tbl}
    fi
    """

# Filter repeatmasker output
# 1) Sequences > 50 nt of consecutive sequence without N
# 2) Sequences with >= 40% of total length of being masked
# input, output, and params names must match function arguments
rule repeatmasker_good:
  input:
    masked = rules.repeatmasker.output.masked,
    original = rules.repeatmasker.input.fa
  output:
    masked_filt = "mask/{run}_repmaskedgood.fa",
    original_filt = "mask/{run}_unmaskedgood.fa"
  params:
    min_length = 50,
    por_n = 40
  wrapper:
    "https://raw.githubusercontent.com/avilab/snakemake-wrappers/master/filter/masked"