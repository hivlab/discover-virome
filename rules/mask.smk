
# Tantan mask of low complexity DNA sequences
rule tantan:
    input:
        rules.fix_fasta.output[0]
    output:
        temp("output/{group}/tantan.fa")
    params:
        extra = "-x N" # mask low complexity using N
    resources:
        runtime = 120,
        mem_mb = 8000
    wrapper:
        f"{WRAPPER_PREFIX}/master/tantan"


# Filter tantan output
# 1) Sequences > 50 nt of consecutive sequence without N
# 2) Sequences with >= 40% of total length of being masked
rule tantan_good:
    input:
        masked = rules.tantan.output[0]
    output:
        masked_filt = temp("output/{group}/tantangood.fa")
    params:
        min_length = 50,
        por_n = 40
    resources:
        runtime = 120
    wrapper:
        f"{WRAPPER_PREFIX}/master/filter/masked"


# Split reads to smaller chunks for RepeatMasker
rule split_fasta:
    input:
        rules.tantan_good.output.masked_filt
    output:
        temp(expand("output/{{group}}/splits/repeatmasker_{n}.fa", n = N))
    params:
        config["splits"]
    resources:
        runtime = lambda wildcards, attempt: 90 + (attempt * 30),
        mem_mb = 4000
    wrapper:
        f"{WRAPPER_PREFIX}/master/split-fasta"


# Repeatmasker
# Outputs are generated from input file names by RepeatMasker
# must have file extension '.masked'
# If no repetitive sequences were detected symlink output to input file
rule repeatmasker:
    input:
        fa = "output/{group}/splits/repeatmasker_{n}.fa"
    output:
        masked = temp("output/{group}/splits/repeatmasker_{n}.fa.masked"),
        out = temp("output/{group}/splits/repeatmasker_{n}.fa.out"),
        cat = temp("output/{group}/splits/repeatmasker_{n}.fa.cat"),
        tbl = "output/{group}/splits/repeatmasker_{n}.fa.tbl"
    params:
        extra = "-qq"
    threads: 8
    resources:
        runtime = 1440,
        mem_mb = 16000
    singularity:
        "shub://tpall/repeatmasker-singularity"
    script:
        f"{WRAPPER_PREFIX}/master/repeatmasker/wrapper.py"


# Filter repeatmasker output
# 1) Sequences > min_length nt of consecutive sequence without N
# 2) Sequences with >= % of total length of being masked
# input, output, and params names must match function arguments
rule repeatmasker_good:
    input:
        masked = rules.repeatmasker.output.masked,
        original = rules.repeatmasker.input.fa
    output:
        masked_filt = temp("output/{group}/splits/repmaskedgood_{n}.fa"),
        original_filt = temp("output/{group}/splits/unmaskedgood_{n}.fa")
    params:
        min_length = 100,
        por_n = 30
    resources:
        runtime = 120
    wrapper:
        f"{WRAPPER_PREFIX}/master/filter/masked"
