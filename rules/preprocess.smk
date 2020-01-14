
FTP = FTPRemoteProvider(username=config["username"], password=config["password"])


def get_fastq(wildcards):
    """Get fraction read file paths from samples.tsv"""
    return list(SAMPLES.loc[(wildcards.group, wildcards.run), ["fq1", "fq2"]])


def get_frac(wildcards):
    """Get fraction of reads to be sampled from samples.tsv"""
    return SAMPLES.loc[(wildcards.group, wildcards.run), ["frac"]][0]


rule preprocess:
    input:
      sample = lambda wildcards: FTP.remote(get_fastq(wildcards), immediate_close=True) if config["remote"] else get_fastq(wildcards)
    output:
      adapters = temp("output/preprocess/{group}-{run}_adapters.fa"),
      merged = temp("output/preprocess/{group}-{run}_merged.fq"),
      unmerged = temp("output/preprocess/{group}-{run}_unmerged.fq"),
      trimmed = temp("output/preprocess/{group}-{run}_trimmed.fq"),
      sampled = temp("output/preprocess/{group}-{run}_sample.fq")
    params:
      bbduk = "qtrim=r trimq=10 maq=10 minlen=100",
      frac = 1, #lambda wildcards: get_frac(wildcards),
      seed = config["seed"]
    threads: 4
    wrapper:
      wrapper_prefix + "master/preprocess"


# Map reads to host.
rule bwa_mem_host:
    input:
      reads = [rules.preprocess.output.sampled]
    output:
      temp("mapped/{group}-{run}_host.bam")
    params:
      db_prefix = HOST_GENOME,
      extra = "-L 100,100 -k 15",
      sorting = "none"
    log:
      "logs/{group}-{run}_bwa_mem_host.log"
    threads: 2
    wrapper:
      "https://raw.githubusercontent.com/tpall/snakemake-wrappers/bug/snakemake_issue145/bio/bwa/mem"


# Extract unmapped reads and convert to fasta.
rule unmapped_host:
    input:
      rules.bwa_mem_host.output
    output:
      fastq = temp("output/preprocess/{group}-{run}_unmapped.fq"),
      fasta = temp("output/preprocess/{group}-{run}_unmapped.fa")
    params:
      reformat_fastq_extra = "-Xmx8000m",
      reformat_fasta_extra = "uniquenames -Xmx8000m"
    wrapper:
      BWA_UNMAPPED


rule assemble:
    input: 
      se = expand("output/preprocess/{{group}}-{run}_unmapped.fq", run = RUN)
    output: 
      contigs = temp("output/{group}/final.contigs.fa")
    params:
      extra = "--min-contig-len 1000"
    threads: 4
    log: "logs/{group}_assemble.log"
    wrapper:
      wrapper_prefix + "release/metformin-pill/assembly/megahit"


# Calculate assembly coverage stats
# nodisk keeps index in memory, otherwise index will be written once to project root (ref/1) from first run to be processed 
# and reused for other unrelated runs.
# Key "input" will be parsed to "in", "input1" to "in1" etc.
rule coverage:
    input:
      ref = rules.assemble.output.contigs, 
      input = expand("output/preprocess/{{group}}-{run}_unmapped.fq", run = RUN) 
    output:
      out = temp("output/contigs/{group}_aln.sam"),
      covstats = "output/stats/{group}_coverage.txt",
      basecov = "output/stats/{group}_basecov.txt"
    params: 
      extra = "kfilter=22 subfilter=15 maxindel=80 nodisk"
    wrapper:
      wrapper_prefix + "master/bbmap/bbwrap"


rule quast:
	input:
		rules.assemble.output.contigs
	output:
		"output/stats/quast/{group}/report.html"
	params:
		outdir = "output/stats/quast/{group}"
	threads: 4
	singularity: 
    "docker://quay.io/biocontainers/quast:5.0.2--py27pl526ha92aebf_0"
	shell: 
    "quast -o {params.outdir} {input} --threads {threads}"


# Run cd-hit to cluster similar contigs
rule cd_hit:
    input:
      rules.assemble.output.contigs
    output:
      repres = temp("output/cdhit/{group}_cdhit.fa"),
      clstr = "output/cdhit/{group}_cdhit.fa.clstr"
    params:
      extra = "-c 0.9 -G 1 -g 1 -prog megablast -s '-num_threads 4'"
    singularity:
      "shub://avilab/singularity-cdhit"
    shell:
      "psi-cd-hit.pl -i {input} -o {output.repres} {params.extra}"


localrules: cleanup
rule cleanup:
    input:
      contigs = rules.assemble.output.contigs,
      repres = rules.cd_hit.output.repres,
      clstr = rules.cd_hit.output.clstr
    output:
      contigs = "output/contigs/{group}_final-contigs.fa",
      repres = "output/contigs/{group}_cdhit.fa",
      clstr = "output/contigs/{group}_cdhit.fa.clstr"
    shell:
      """
      mv {input.contigs} {output.contigs}
      mv {input.repres} {output.repres}
      mv {input.clstr} {output.clstr}
      rm -rf $(dirname {input})
      """


# Tantan mask of low complexity DNA sequences
rule tantan:
    input:
      rules.cleanup.output.repres
    output:
      temp("output/RM/{group}_tantan.fasta")
    params:
      extra = "-x N" # mask low complexity using N
    wrapper:
      wrapper_prefix + "master/tantan"


# Filter tantan output
# 1) Sequences > 50 nt of consecutive sequence without N
# 2) Sequences with >= 40% of total length of being masked
rule tantan_good:
    input:
      masked = rules.tantan.output
    output:
      masked_filt = temp("output/RM/{group}_repeatmasker.fa")
    params:
      min_length = 50,
      por_n = 40
    wrapper:
      LN_FILTER


# Repeatmasker
# Outputs are generated from input file names by RepeatMasker
# must have file extension '.masked'
# If no repetitive sequences were detected symlink output to input file
rule repeatmasker:
    input:
      fa = rules.tantan_good.output
    output:
      masked = temp("output/RM/{group}_repeatmasker.fa.masked"),
      out = temp("output/RM/{group}_repeatmasker.fa.out"),
      cat = temp("output/RM/{group}_repeatmasker.fa.cat"),
      tbl = "output/RM/{group}_repeatmasker.fa.tbl"
    params:
      extra = "-qq"
    threads: 8
    singularity:
      "shub://tpall/repeatmasker-singularity"
    script:
      RM


# Filter repeatmasker output
# 1) Sequences > 50 nt of consecutive sequence without N
# 2) Sequences with >= 40% of total length of being masked
# input, output, and params names must match function arguments
rule repeatmasker_good:
    input:
      masked = rules.repeatmasker.output.masked,
      original = rules.tantan_good.output
    output:
      masked_filt = temp("output/RM/{group}_repmaskedgood.fa"),
      original_filt = temp("output/RM/{group}_unmaskedgood.fa")
    params:
      min_length = 50,
      por_n = 40
    wrapper:
      LN_FILTER


# Split reads to smaller chunks
rule split_fasta:
    input:
      rules.repeatmasker_good.output.masked_filt
    output:
      temp(expand("output/RM/{{group}}_repmaskedgood_{n}.fa", n = N))
    params:
      config["split_fasta"]["n_files"]
    wrapper:
      wrapper_prefix + "master/split-fasta"


# Collect stats from preprocess outputs.
rule preprocess_stats:
    input:
      expand("output/preprocess/{{group}}-{run}_{file}.fq", run = RUN, file = ["trimmed", "unmapped"]),
      rules.cleanup.output.contigs,
      rules.cd_hit.output.repres,
      rules.tantan.output,
      rules.tantan_good.output,
      rules.repeatmasker_good.output
    output:
      "output/stats/{group}_preprocess-stats.tsv"
    params:
      extra = "-T"
    wrapper:
      SEQ_STATS


# host mapping stats.
rule host_bam_stats:
    input:
      rules.bwa_mem_host.output
    output:
      "output/stats/{group}-{run}_host-bam-stats.txt"
    params:
      extra = "-f 4",
      region = ""
    wrapper:
        "0.42.0/bio/samtools/stats"
