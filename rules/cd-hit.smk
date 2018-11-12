
## Run cd-hit to find and munge duplicate reads
rule cd_hit:
  input: rules.refgenomefilter.output.fq
  output:
    fa = "avasta/refgenomefilter/{sample}_refgenome_unmapped.fa",
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
    sed -n '1~4s/^@/>/p;2~4p' {input} > {output.fa}
    cd-hit-est -i {output.fa} -o {output.repres} -T {threads} {params} > {output.report}
    """

## Add top n (3) similar sequences to cluster representative sequences
## outputs file with sequence ids (clstr) to subset refgenome unmapped output
## fasta file and appends these sequences to cluster representative sequences (repres)
rule parse_cdhit:
  input:
      repres = rules.cd_hit.output.repres,
      clstr = rules.cd_hit.output.clstr,
      fq = rules.refgenomefilter.output.fq,
      join = rules.fastq_join.output[0]
  output:
      topn_clstr = "avasta/cdhit/{sample}_cdhit_topn.clstr",
      topn = "avasta/cdhit/{sample}_cdhit_topn.fa",
      merged = "avasta/cdhit/{sample}_cdhit_merged.fa",
      join = "avasta/cdhit/{sample}_cdhit_join.fa",
      un = "avasta/cdhit/{sample}_cdhit_un.fa"
  params:
      top_n = 3
  conda:
      "../envs/scipy.yml"
  script:
      "../scripts/parse_cdhit.py"
