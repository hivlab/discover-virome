# Metagenome assembly
rule assemble:
    input: 
      pe1 = rules.refgenome_unmapped_bam2fq.output[0],
      pe2 = rules.refgenome_unmapped_bam2fq.output[1]
    output: 
      contigs = "assemble/{run}/final.contigs.fa"
    params:
      options = "--min-contig-len 500"
    threads: 2
    log: "logs/{run}_assemble.log"
    wrapper:
      "https://bitbucket.org/tpall/snakemake-wrappers/raw/77183b4bdef5103a2c3e60d4d6c3825a17d5debc/bio/assembly/megahit/"

# Calculate assembly coverage stats
rule coverage:
    input: 
      ref = rules.assemble.output.contigs,
      in1 = rules.refgenome_unmapped_bam2fq.output[0],
      in2 = rules.refgenome_unmapped_bam2fq.output[1]
    output:
      aln = temp("assemble/{run}/aln.sam.gz"),
      cov = "assemble/{run}/coverage.txt"
    params: 
      options = "kfilter=22 subfilter=15 maxindel=80"
    wrapper:
      "https://bitbucket.org/tpall/snakemake-wrappers/raw/99269f6f66af772fd045a05435915d6d7c9c3533/bio/assembly/coverage"
