
## Run cd-hit to find and munge duplicate reads
rule cd_hit:
  input: rules.refgenome_unmapped.output.fa
  output:
    repres = "cdhit/{sample}_cdhit.fa",
    clstr = "cdhit/{sample}_cdhit.fa.clstr",
    report = "cdhit/{sample}_cdhit.report"
  params:
    "-c 0.95 -G 0 -n 8 -d 0 -aS 0.95 -g 1 -r 1 -M 0"
  threads: 8
  conda:
    "../envs/cd-hit.yml"
  shell:
    """
    cd-hit-est -i {input} -o {output.repres} -T {threads} {params} > {output.report}
    """
## Get top n (3) similar sequences to cluster representative
rule parse_cd_hit:
  input:
      rules.cd_hit.output.clstr,
      rules.refgenome_unmapped.output.fa
  output:
      "cdhit/{sample}_cdhit_clstr.topn",
      "cdhit/{sample}_cdhit_topn.fa"
  conda:
      "../envs/bioconda.yml"
  script:
      "../scripts/parse_clstr.py"
