# Runs spades metagenome+assembly module
rule assemble:
    input: 
      rules.fastp.output
    output: 
      contigs = "assemble/{sample}.contigs.fa"
    params:
      options = "--out-prefix {sample} --min-contig-len 500"
    threads: 8
    log: "logs/{sample}_assemble.log"
    wrapper:
      "https://raw.githubusercontent.com/avilab/snakemake-wrappers/master/assembly/megahit"