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
rule bwa_mem:
  input:
    reads = [rules.fastp.output]
  output:
    "mapped/{run}_refgenome_mapped.bam"
  params:
    index = config["ref_genome"],
    sort = "samtools",
    sort_order = "queryname"
  log:
    "logs/{run}_bwa_mem_refgenome.log"
  threads: 2
  wrapper:
    "0.32.0/bio/bwa/mem"

rule samtools_view:
  input:
    rules.bwa_mem.output
  output:
    "mapped/{run}_refgenome_unmapped.bam"
  params:
    "-b -f 4" # bam output
  wrapper:
    "0.32.0/bio/samtools/view"

# Convert bam file to fastq files.
rule bamtofastq:
  input:
    rules.samtools_view.output
  output:
    "mapped/{run}_unmapped_1.fq",
    "mapped/{run}_unmapped_2.fq"
  log: 
    "logs/{run}_bamtofastq.log"
  wrapper:
    "https://bitbucket.org/tpall/snakemake-wrappers/raw/8e23fd260cdbed02450a7eb1796dce984d2e1f8f/bio/bedtools/bamtofastq"
