FTP = FTPRemoteProvider(username = config["username"], password = config["password"])

def get_fastq(wildcards):
    """Get read file urls from samples.tsv"""
    urls = SAMPLES.loc[wildcards.run, ['fq1', 'fq2']]
    return list(urls)

# Adapter trimming and quality filtering.
rule fastp:
  input:
    lambda wildcards: FTP.remote(get_fastq(wildcards), immediate_close = True) if config["remote"] else get_fastq(wildcards)
  output:
    "trim/{run}_read1_trimmed.fq.gz",
    "trim/{run}_read2_trimmed.fq.gz"
  params:
    options = "--trim_front1 5 --trim_tail1 5 --length_required 50 --low_complexity_filter --complexity_threshold 8",
    html = temp("munge/{run}_fastp_report.html"),
    json = "munge/{run}_fastp_report.json"
  threads: 2
  wrapper:
    "https://bitbucket.org/tpall/snakemake-wrappers/raw/8e23fd260cdbed02450a7eb1796dce984d2e1f8f/bio/fastp"

# Align sequences to reference genome and extract unmapped reads.
rule refgenome_unmapped:
  input:
    index = config["ref_genome"],
    reads = rules.fastp.output
  output:
    "trim/{run}_refgenome_unmapped.bam"
  params:
    sort = "unmapped",
    sort_order = "queryname"
  log:
    "logs/{run}_bwa_map_refgenome.log"
  threads: 2
  wrapper:
    "https://bitbucket.org/tpall/snakemake-wrappers/raw/8e23fd260cdbed02450a7eb1796dce984d2e1f8f/bio/bwa/mem"

# Convert bam file to fastq files.
rule bamtofastq:
  input:
    rules.refgenome_unmapped.output
  output:
    "trim/{run}_unmapped_1.fq",
    "trim/{run}_unmapped_2.fq"
  log: 
    "logs/{run}_bamtofastq.log"
  wrapper:
    "https://bitbucket.org/tpall/snakemake-wrappers/raw/8e23fd260cdbed02450a7eb1796dce984d2e1f8f/bio/bedtools/bamtofastq"
