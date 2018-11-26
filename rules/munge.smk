
def get_fastq(wildcards):
 return SAMPLES.loc[wildcards.sample, ['fq1', 'fq2']].dropna()

## Preprocessing for fastq files
# Adapter trimming
# Quality filtering
# Replaces AdapteRemoval, prinseq and fastqc
rule fastp:
    input:
      get_fastq
    output:
      pair1 = "munge/{sample}_pair1_trimmed.gz",
      pair2 = "munge/{sample}_pair2_trimmed.gz",
      html = "munge/{sample}_fastp_report.html",
      json = "munge/{sample}_fastp_report.json"
    params:
      "-f 5 -t 5 -l 50 -y -Y 8"
    threads: 8
    conda:
      "../envs/fastp.yaml"
    log: "logs/{sample}_fastp.log"
    shell:
      """
      fastp -i {input[0]} -I {input[1]} -o {output.pair1} -O {output.pair2} {params} -h {output.html} -j {output.json} -w {threads} > {log} 2>&1
      """

## Stitch paired reads
rule fastq_join:
  input:
    rules.fastp.output.pair1,
    rules.fastp.output.pair2
  output:
    "munge/{sample}_join.fq.gz",
    "munge/{sample}_un1.fq.gz",
    "munge/{sample}_un2.fq.gz"
  params:
    diff = config["fastq-join"]["maximum_difference"],
    overlap = config["fastq-join"]["minimum_overlap"],
    template = "munge/{sample}_%.fq.gz"
  conda:
    "../envs/fastq-join.yaml"
  log:
    "logs/{sample}_fastq_join.log"
  shell:
    """
    fastq-join -p {params.diff} -m {params.overlap} {input} -o {params.template} > {log} 2>&1
    """
