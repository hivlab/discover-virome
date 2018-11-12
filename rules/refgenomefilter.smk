
## Align sequences to reference genome and extract unmapped reads
rule refgenome_unmapped:
    input:
        db = config["ref_genome"],
        reads = rules.fastq_join.output
    output:
      merged = "munge/{sample}_merge_reads.fq.gz",
      bam = "avasta/refgenomefilter/{sample}_refgenome_unmapped.bam",
      fq = "avasta/refgenomefilter/{sample}_refgenome_unmapped.fq"
    log:
        "logs/avasta/{sample}_bwa_map_refgenome.log"
    threads: 8
    conda:
      "../envs/bwa-sam-bed.yml"
    shell:
      """
      cat {input.reads} > {output.merged}
      (bwa mem -L 100,100 -k 15 -t {threads} {input.db} {ouput.merged} | samtools view -b -S -f 4 - > {output.bam}) 2> {log}
      bedtools bamtofastq -i {output.bam} -fq {output.fq}
      """
