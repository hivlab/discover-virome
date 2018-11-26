## Run cd-hit to define unique sequences
rule cd_hit:
  input:
    rules.filter_contigs.output
  output:
    repres = "avasta/cdhit/{sample}_cdhit.fa",
    report = "avasta/cdhit/{sample}_cdhit.report",
    clstr = "avasta/cdhit/{sample}_cdhit.fa.clstr"
  params:
    "-c 0.95 -G 0 -n 8 -d 0 -aS 0.95 -g 1 -r 1 -M 0"
  threads: 8
  conda:
    "../envs/cd-hit.yaml"
  shell:
    """
    cd-hit-est -i {input} -o {output.repres} -T {threads} {params} > {output.report}
    """
