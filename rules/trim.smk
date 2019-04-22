FTP = FTPRemoteProvider(username = config["username"], password = config["password"])

def get_fastq(wildcards):
    """Get fraction read file paths from samples.tsv"""
    urls = SAMPLES.loc[wildcards.sample, ['fq1', 'fq2']]
    return list(urls)

def get_frac(wildcards):
    """Get fraction of reads to be sampled from samples.tsv"""
    frac = SAMPLES.loc[wildcards.sample, ['frac']][0]
    return frac

# Adapter trimming and quality filtering.
rule fastp:
  input:
    lambda wildcards: FTP.remote(get_fastq(wildcards), immediate_close=True) if config["remote"] else get_fastq(wildcards)
  output:
    "munge/{sample}_read1_trimmed.fq.gz",
    "munge/{sample}_read2_trimmed.fq.gz"
  params:
    options = "--trim_front1 5 --trim_tail1 5 --length_required 50 --low_complexity_filter --complexity_threshold 8",
    html = temp("munge/{sample}_fastp_report.html"),
    json = "munge/{sample}_fastp_report.json"
  threads: 8
  wrapper:
    "https://bitbucket.org/tpall/snakemake-wrappers/raw/8e23fd260cdbed02450a7eb1796dce984d2e1f8f/bio/fastp"

# Stitch paired reads.
rule fastq_join:
  input:
    rules.fastp.output
  output:
    "munge/{sample}_un1.fq.gz",
    "munge/{sample}_un2.fq.gz",
    "munge/{sample}_join.fq.gz"
  params:
    options = "-p 5 -m 10"
  log: "logs/{sample}_fastq_join.log"
  wrapper:
    "https://bitbucket.org/tpall/snakemake-wrappers/raw/8e23fd260cdbed02450a7eb1796dce984d2e1f8f/bio/fastq-join"

## Align sequences to reference genome and extract unmapped reads
rule refgenome_unmapped:
    input:
        config["ref_genome"],
        [rules.fastq_join.output]
    output:
      bam = "refgenomefilter/{sample}_refgenome_unmapped_{n}.bam",
      fq = "refgenomefilter/{sample}_refgenome_unmapped_{n}.fq",
      fa = "refgenomefilter/{sample}_refgenome_unmapped_{n}.fa"
    log:
        "logs/{sample}_bwa_map_refgenome_{n}.log"
    threads: 8
    conda:
      "../envs/bwa-sam-bed.yaml"
    shell:
      """
        (bwa mem -L 100,100 -k 15 -t {threads} {input} | samtools view -b -S -f 4 - > {output.bam}) 2> {log}
        bedtools bamtofastq -i {output.bam} -fq {output.fq}
        cat {output.fq} | sed -n '1~4s/^@/>/p;2~4p' > {output.fa}
      """