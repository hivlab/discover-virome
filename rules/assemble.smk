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
