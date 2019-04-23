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
    "trim/{run}_trimmed_1.fq.gz",
    "trim/{run}_trimmed_2.fq.gz"
  params:
    options = "--trim_front1 5 --trim_tail1 5 --length_required 50 --low_complexity_filter --complexity_threshold 8",
    html = temp("trim/{run}_fastp_report.html"),
    json = "trim/{run}_fastp_report.json"
  threads: 2
  wrapper:
    "https://bitbucket.org/tpall/snakemake-wrappers/raw/8e23fd260cdbed02450a7eb1796dce984d2e1f8f/bio/fastp"

# Align sequences to reference genome and extract unmapped reads.
rule bwa_mem_refgenome:
  input:
    reads = [rules.fastp.output]
  output:
    "trim/{run}_refgenome_mapped.bam"
  params:
    index = config["ref_genome"],
    sort = "none"
  log:
    "logs/{run}_bwa_mem_refgenome.log"
  threads: 2
  wrapper:
    "0.32.0/bio/bwa/mem"

# Convert bam to fastq files
rule refgenome_unmapped_bam2fq:
    input:
      rules.bwa_mem_refgenome.output
    output:
      "trim/{run}_unmapped_1.fq",
      "trim/{run}_unmapped_2.fq"
    params:
      sort = "",
      bam2fq = "-n -f 4"
    threads: 2
    wrapper:
      "0.32.0/bio/samtools/bam2fq/separate"

# Calculate bam stats
rule refgenome_stats:
    input:
      "trim/{run}_refgenome_mapped.bam"
    output:
      "trim/{run}_refgenome_stats.txt"
    params:
      extra = "-f 4",
      region = ""
    log:
      "logs/{run}_refgenome_stats.log"
    wrapper:
        "0.32.0/bio/samtools/stats"
