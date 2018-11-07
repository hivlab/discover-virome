
## Align sequences to reference genome and extract unmapped reads
rule refgenome_unmapped:
    input:
        config["ref_genome"],
        [rules.merge_reads.output.fa]
    output:
      bam = "avasta/refgenomefilter/{sample}_refgenome_unmapped.bam",
      fq = "avasta/refgenomefilter/{sample}_refgenome_unmapped.fq",
      fa = "avasta/refgenomefilter/{sample}_refgenome_unmapped.fa"
    log:
        "logs/avasta/{sample}_bwa_map_refgenome.log"
    threads: 8
    conda:
      "../envs/bwa-sam-bed.yml"
    shell:
      """
        (bwa mem -L 100,100 -k 15 -t {threads} {input} | samtools view -b -S -f 4 - > {output.bam}) 2> {log}
        bedtools bamtofastq -i {output.bam} -fq {output.fq}
        cat {output.fq} | sed -n '1~4s/^@/>/p;2~4p' > {output.fa}
      """
