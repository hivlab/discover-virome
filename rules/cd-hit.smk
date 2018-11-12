
## Run cd-hit to find and munge duplicate reads
rule cd_hit:
  input: rules.refgenomefilter.output.fq
  output:
    fa = "munge/{sample}_merge_reads.fasta",
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
    zcat {input} | sed -n '1~4s/^@/>/p;2~4p' > {output.fa}
    cd-hit-est -i {output.fa} -o {output.repres} -T {threads} {params} > {output.report}
    """

## Add top n (3) similar sequences to cluster representative sequences
## outputs file with sequence ids (clstr) to subset refgenome unmapped output
## fasta file and appends these sequences to cluster representative sequences (repres)
rule parse_cdhit:
  input:
      repres = rules.cd_hit.output.repres,
      clstr = rules.cd_hit.output.clstr,
      fq = rules.refgenomefilter.output.fq
  output:
      topn_clstr = "avasta/cdhit/{sample}_cdhit_topn.clstr",
      topn_fq = "avasta/cdhit/{sample}_cdhit_topn.fq"
  params:
      top_n = 3
  conda:
      "../envs/scipy.yml"
  script:
      "../scripts/parse_cdhit.py"
