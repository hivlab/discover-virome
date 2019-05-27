
FTP = FTPRemoteProvider(username = config["username"], password = config["password"])

def get_fastq(wildcards):
    """Get fraction read file paths from samples.tsv"""
    urls = SAMPLES.loc[wildcards.run, ['fq1', 'fq2']]
    return list(urls)

def get_frac(wildcards):
    """Get fraction of reads to be sampled from samples.tsv"""
    frac = SAMPLES.loc[wildcards.run, ['frac']][0]
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
    reformat_fastq_extra = "-Xmx8000m",
    reformat_fasta_extra = "uniquenames"
  wrapper:
    "https://raw.githubusercontent.com/avilab/vs-wrappers/master/unmapped"

rule assemble:
    input: 
      pe12 = rules.unmapped_refgenome.output.fastq,
    output: 
      contigs = "assemble/{run}/final.contigs.fa"
    params:
      options = "--min-contig-len 500"
    threads: 2
    log: "logs/{run}_assemble.log"
    wrapper:
      "https://bitbucket.org/tpall/snakemake-wrappers/raw/77183b4bdef5103a2c3e60d4d6c3825a17d5debc/bio/assembly/megahit"

# Calculate assembly coverage stats
rule bbwrap:
    input:
      ref = rules.assemble.output.contigs, 
      input = rules.unmapped_refgenome.output.fastq # input will be parsed to 'in', input1 to in1 etc.
    output:
      out = pipe("assemble/{run}/aln.sam")
    wrapper:
      "https://raw.githubusercontent.com/avilab/vs-wrappers/master/bbmap/bbwrap"

rule coverage:
    input: 
      input = "assemble/{run}/aln.sam" # input will be parsed to 'in', input1 to in1 etc.
    output:
      cov = "assemble/stats/{run}_coverage.txt"
    params: 
      extra = "kfilter=22 subfilter=15 maxindel=80"
    wrapper:
      "https://raw.githubusercontent.com/avilab/vs-wrappers/master/bbmap/pileup"

# Tantan mask of low complexity DNA sequences
rule tantan:
  input:
    rules.assemble.output
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
    out = temp("assemble/mask/{run}_repeatmasker.fa.out"),
    cat = temp("assemble/mask/{run}_repeatmasker.fa.cat"),
    tbl = "assemble/mask/{run}_repeatmasker.fa.tbl"
  params:
    outdir = "assemble/mask"
  threads: 2
  singularity:
    "shub://tpall/repeatmasker-singularity"
  shell:
    """
    RepeatMasker -qq -pa {threads} {input} -dir {params.outdir}
    if head -n 1 {output.out} | grep -q "There were no repetitive sequences detected"
      then ln -sr {input} {output.masked} \
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
    original = rules.tantan_good.output
  output:
    masked_filt = temp("assemble/mask/{run}_repmaskedgood.fa"),
    original_filt = temp("assemble/mask/{run}_unmaskedgood.fa")
  params:
    min_length = 50,
    por_n = 40
  wrapper:
    "https://raw.githubusercontent.com/avilab/snakemake-wrappers/master/filter/masked"

# MegaBlast against reference genome to remove host sequences
rule megablast_refgenome:
    input:
      query = rules.repeatmasker_good.output.masked_filt
    output:
      out = temp("assemble/blast/{run}_megablast.tsv")
    params:
      db = config["ref_genome"],
      task = "megablast",
      perc_identity = config["megablast_ref_genome"]["perc_identity"],
      evalue = config["megablast_ref_genome"]["evalue"],
      word_size = config["megablast_ref_genome"]["word_size"],
      max_hsps = config["blastn_virus"]["max_hsps"],
      show_gis = True,
      num_threads = 2,
      outfmt = "'6 qseqid sgi pident length mismatch gapopen qstart qend sstart send evalue bitscore'"
    wrapper:
      config["wrappers"]["blast"]

# Filter megablast records for the cutoff value
rule parse_megablast:
    input:
      blast_result = rules.megablast_refgenome.output.out,
      query = rules.repeatmasker_good.output.masked_filt
    output:
      mapped = temp("assemble/blast/{run}_refgenome_megablast_known-host.tsv"),
      unmapped = temp("assemble/blast/{run}_refgenome_megablast_unmapped.fa")
    params:
      e_cutoff = 1e-10,
      outfmt = rules.megablast_refgenome.params.outfmt
    wrapper:
      config["wrappers"]["parse_blast"]

# Collect stats from preprocess outputs.
rule preprocess_stats:
  input:
    rules.preprocess.output.trimmed,
    rules.unmapped_refgenome.output,
    rules.parse_megablast.output.unmapped,
    rules.tantan.output,
    rules.tantan_good.output,
    rules.repeatmasker_good.output
  output:
    "assemble/stats/{run}_preprocess.tsv"
  params:
    extra = "-T"
  wrapper:
    config["wrappers"]["stats"]

# Refgenome mapping stats.
rule refgenome_bam_stats:
    input:
      rules.bwa_mem_refgenome.output
    output:
      "assemble/stats/{run}_refgenome_stats.txt"
    params:
      extra = "-f 4",
      region = ""
    wrapper:
        "0.32.0/bio/samtools/stats"
