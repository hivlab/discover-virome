
## Run cd-hit to find and munge duplicate reads
rule cd_hit:
  input: rules.refgenome_unmapped.output.fa
  output:
    repres = "avasta/cdhit/{sample}_cdhit.fa",
    clstr = "avasta/cdhit/{sample}_cdhit.fa.clstr",
    report = "avasta/cdhit/{sample}_cdhit.report"
  params:
    "-c 0.95 -G 0 -n 8 -d 0 -aS 0.95 -g 1 -r 1 -M 0"
  threads: 8
  conda:
    "../envs/cd-hit.yml"
  shell:
    """
    cd-hit-est -i {input} -o {output.repres} -T {threads} {params} > {output.report}
    """
## Add top n (3) similar sequences to cluster representative sequences
## outputs file with sequence ids (clstr) to subset refgenome unmapped output
## fasta file and appends these sequences to cluster representative sequences (repres)
rule parse_cdhit:
  input:
      repres = rules.cd_hit.output.repres,
      clstr = rules.cd_hit.output.clstr,
      fa = rules.refgenome_unmapped.output.fa
  output:
      topn_clstr = "avasta/cdhit/{sample}_cdhit_topn.clstr",
      topn_fa = "avasta/cdhit/{sample}_cdhit_topn.fa"
  params:
      top_n = 3
  conda:
      "../envs/bioconda.yml"
  script:
      "../scripts/parse_cdhit.py"
