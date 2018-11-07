
## Align sequences to reference genome and extract unmapped reads
rule refgenome_unmapped:
    input:
        config["ref_genome"],
        [rules.merge_reads.output.fa]
    output:
      bam = "refgenomefilter/{sample}_refgenome_unmapped_{n}.bam",
      fq = "refgenomefilter/{sample}_refgenome_unmapped_{n}.fq",
      fa = "refgenomefilter/{sample}_refgenome_unmapped_{n}.fa"
    log:
        "logs/{sample}_bwa_map_refgenome_{n}.log"
    threads: 8
    conda:
      "../envs/bwa-sam-bed.yml"
    shell:
      """
        (bwa mem -L 100,100 -k 15 -t {threads} {input} | samtools view -b -S -f 4 - > {output.bam}) 2> {log}
        bedtools bamtofastq -i {output.bam} -fq {output.fq}
        cat {output.fq} | sed -n '1~4s/^@/>/p;2~4p' > {output.fa}
      """
