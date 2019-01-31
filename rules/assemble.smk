# Runs spades metagenome+assembly module
rule assemble:
    input: 
      pe1 = rules.fastp.output[0],
      pe2 = rules.fastp.output[1]
    output: 
      contigs = "assemble/{sample}/final.contigs.fa"
    params:
      options = "--min-contig-len 500"
    threads: 8
    log: "logs/{sample}_assemble.log"
    wrapper:
      "https://raw.githubusercontent.com/avilab/snakemake-wrappers/master/assembly/megahit"

rule align:
    input: 
      ref = rules.assemble.output.contigs,
      in1 = rules.fastp.output[0],
      in2 = rules.fastp.output[1]
    output:
      out = "align/{sample}/aln.sam.gz",
      path = directory("index/{sample}")
    params:
      options = "kfilter=22 subfilter=15 maxindel=80"
    wrapper:
      "https://raw.githubusercontent.com/avilab/snakemake-wrappers/master/bbmap/bbwrap"

rule coverage:
    input: 
      rules.align.output.out
    output:
      out = "align/{sample}/coverage.tsv"
    log: "logs/{sample}_coverage.log"
    wrapper:
      "https://raw.githubusercontent.com/avilab/snakemake-wrappers/master/bbmap/pileup"

rule samtools_view:
    input:
        rules.align.output.out
    output:
        "mapped/{sample}.bam"
    params:
        "-b -F 4 -q 10"
    wrapper:
        "0.31.1/bio/samtools/view"

rule samtools_sort:
    input:
        rules.samtools_view.output
    output:
        "mapped/sorted_{sample}.bam"
    params:
        "-n"
    threads: 8
    wrapper:
        "0.31.1/bio/samtools/sort"

