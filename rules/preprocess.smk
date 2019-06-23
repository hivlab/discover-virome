
FTP = FTPRemoteProvider(username = config["username"], password = config["password"])

def get_fastq(wildcards):
    """Get fraction read file paths from samples.tsv"""
    urls = RUNS.loc[wildcards.run, ['fq1', 'fq2']]
    return list(urls)

def get_frac(wildcards):
    """Get fraction of reads to be sampled from samples.tsv"""
    frac = RUNS.loc[wildcards.run, ['frac']][0]
    return frac

rule preprocess:
  input:
    sample = lambda wildcards: FTP.remote(get_fastq(wildcards), immediate_close=True) if config["remote"] else get_fastq(wildcards)
  output:
    adapters = temp("preprocess/{run}_adapters.fa"),
    merged = temp("preprocess/{run}_merged.fq"),
    unmerged = temp("preprocess/{run}_unmerged.fq"),
    reads = temp("preprocess/{run}_reads.fq"),
    trimmed = temp("preprocess/{run}_trimmed.fq"),
    sampled = temp("preprocess/{run}_sample.fq")
  params:
    bbduk = "qtrim=r trimq=10 maq=10 minlen=100",
    frac = lambda wildcards: get_frac(wildcards),
    seed = config["seed"]
  threads: 2
  wrapper:
    "https://raw.githubusercontent.com/avilab/vs-wrappers/master/preprocess"

# Map reads to Refgenome.
rule bwa_mem_refgenome:
  input:
    reads = [rules.preprocess.output.sampled]
  output:
    temp("mapped/{run}_refgenome.bam")
  params:
    index = config["ref_genome"],
    extra = "-L 100,100 -k 15",
    sort = "none"
  log:
    "logs/{run}_bwa_map_refgenome.log"
  threads: 2
  wrapper:
    "0.32.0/bio/bwa/mem"

# Extract unmapped reads and convert to fasta.
rule unmapped_refgenome:
  input:
    rules.bwa_mem_refgenome.output
  output:
    fastq = temp("preprocess/{run}_unmapped.fq"),
    fasta = temp("preprocess/{run}_unmapped.fa")
  params:
    reformat_fastq_extra = "-Xmx2000m",
    reformat_fasta_extra = "uniquenames -Xmx2000m"
  wrapper:
    "https://raw.githubusercontent.com/avilab/vs-wrappers/master/unmapped"

rule assemble:
  input: 
    se = rules.unmapped_refgenome.output.fastq,
  output: 
    contigs = temp("assemble/{run}/final.contigs.fa"),
    intermediate_contigs = temp(directory("assemble/{run}"))
  params:
    options = "--min-contig-len 500"
  threads: 2
  log: "logs/{run}_assemble.log"
  wrapper:
    "https://bitbucket.org/tpall/snakemake-wrappers/raw/adc9201669a4c121968ac044ad149e9b292774d8/bio/assembly/megahit"

# Calculate assembly coverage stats
# nodisk keeps index in memory, otherwise index will be written once to project root (ref/1) from first run to be processed 
# and reused for other unrelated runs
rule bbwrap:
  input:
    ref = rules.assemble.output.contigs, 
    input = rules.unmapped_refgenome.output.fastq # input will be parsed to 'in', input1 to in1 etc.
  output:
    out = temp("assemble/{run}/aln.sam")
  params: 
    extra = "kfilter=22 subfilter=15 maxindel=80 nodisk"
  wrapper:
    "https://raw.githubusercontent.com/avilab/vs-wrappers/master/bbmap/bbwrap"

rule coverage:
  input: 
    input = rules.bbwrap.output # input will be parsed to 'in', input1 to in1 etc.
  output:
    out = "assemble/stats/{run}_assembly-coverage.txt"
  wrapper:
    "https://raw.githubusercontent.com/avilab/vs-wrappers/master/bbmap/pileup"

# Filter contigs by setting minimum threshold for average coverage
rule coverage_good:
  input:
    contigs = rules.assemble.output.contigs,
    coverage = rules.coverage.output.out
  output:
    temp("assemble/{run}/good_contigs.fa")
  params:
    avg_coverage = 8 # average coverage threshold 
  wrapper:
    "https://raw.githubusercontent.com/avilab/vs-wrappers/master/assembly/filter_coverage"

# Tantan mask of low complexity DNA sequences
rule tantan:
  input:
    rules.coverage_good.output
  output:
    temp("assemble/mask/{run}_tantan.fasta")
  params:
    extra = "-x N" # mask low complexity using N
  wrapper:
    "https://bitbucket.org/tpall/snakemake-wrappers/raw/7e681180a5607f20594b3070f8eced7ccd245a89/bio/tantan"

# Filter tantan output
# 1) Sequences > 50 nt of consecutive sequence without N
# 2) Sequences with >= 40% of total length of being masked
rule tantan_good:
  input:
    masked = rules.tantan.output
  output:
    masked_filt = temp("assemble/mask/{run}_repeatmasker.fa")
  params:
    min_length = 50,
    por_n = 40
  wrapper:
    "https://raw.githubusercontent.com/avilab/snakemake-wrappers/master/filter/masked"

# Repeatmasker
# Outputs are generated from input file names by RepeatMasker
# must have file extension '.masked'
# If no repetitive sequences were detected symlink output to input file
rule repeatmasker:
  input:
    rules.tantan_good.output
  output:
    masked = temp("assemble/mask/{run}_repeatmasker.fa.masked"),
    out = temp("assemble/mask/{run}_repeatmasker.fa.out")
  shadow: "full"
  params:
    outdir = "assemble/mask"
  threads: 2
  singularity:
    "shub://tpall/repeatmasker-singularity"
  shell:
    """
    RepeatMasker -qq -pa {threads} {input} -dir {params.outdir}
    if head -n 1 {output.out} | grep -q "There were no repetitive sequences detected"
      then ln -sr {input} {output.masked}
    fi
    """

# Filter repeatmasker output
# 1) Sequences > 50 nt of consecutive sequence without N
# 2) Sequences with >= 40% of total length of being masked
# input, output, and params names must match function arguments
rule repeatmasker_good:
  input:
    masked = rules.repeatmasker.output.masked,
    original = rules.tantan_good.output
  output:
    masked_filt = temp("assemble/mask/{run}_repmaskedgood.fa"),
    original_filt = temp("assemble/mask/{run}_unmaskedgood.fa")
  params:
    min_length = 50,
    por_n = 40
  wrapper:
    "https://raw.githubusercontent.com/avilab/snakemake-wrappers/master/filter/masked"

# Split reads to smaller chunks for Repeatmasker
rule split_fasta:
  input:
    rules.repeatmasker_good.output.masked_filt
  output:
    temp(expand("assemble/mask/{{run}}_repmaskedgood_{n}.fa", n = N))
  params:
    config["split_fasta"]["n_files"]
  wrapper:
    "https://bitbucket.org/tpall/snakemake-wrappers/raw/7e681180a5607f20594b3070f8eced7ccd245a89/bio/split-fasta"

# Collect stats from preprocess outputs.
rule preprocess_stats:
  input:
    rules.preprocess.output.trimmed,
    rules.assemble.output.contigs,
    rules.coverage_good.output,
    rules.unmapped_refgenome.output,
    rules.tantan.output,
    rules.tantan_good.output,
    rules.repeatmasker_good.output
  output:
    "assemble/stats/{run}_assembly-preprocess.tsv"
  params:
    extra = "-T"
  wrapper:
    config["wrappers"]["stats"]

# Refgenome mapping stats.
rule refgenome_bam_stats:
    input:
      rules.bwa_mem_refgenome.output
    output:
      "assemble/stats/{run}_assembly-refgenome-stats.txt"
    params:
      extra = "-f 4",
      region = ""
    wrapper:
        "0.32.0/bio/samtools/stats"
