
FTP = FTPRemoteProvider(username = config["username"], password = config["password"])

def get_fastq(wildcards):
    """Get fraction read file paths from samples.tsv"""
    urls = RUNS.loc[wildcards.run, ['fq1', 'fq2']]
    return list(urls)

def get_frac(wildcards):
    """Get fraction of reads to be sampled from samples.tsv"""
    frac = RUNS.loc[wildcards.run, ['frac']][0]
    return frac

# Convert reads to interleaved format
rule interleave:
    input:
        lambda wildcards: FTP.remote(get_fastq(wildcards), immediate_close=True) if config["remote"] else get_fastq(wildcards)
    output:
        out = temp("output/{run}/interleaved.fq.gz"),
        bhist = "output/{run}/bhist.txt",
        qhist = "output/{run}/qhist.txt",
        aqhist = "output/{run}/aqhist.txt",
        bqhist = "output/{run}/bqhist.txt",
        lhist = "output/{run}/lhist.txt",
        gchist = "output/{run}/gchist.txt"
    params:
        extra = "-Xmx4g"
    resources:
        runtime = 120,
        mem_mb = 4000
    log:
        "output/{run}/log/interleave.txt"
    wrapper:
        WRAPPER_PREFIX + "master/bbtools/reformat"


# Remove PCR and optical duplicates
rule clumpify:
    input:
        input = rules.interleave.output.out
    output:
        out = temp("output/{run}/clumpify.fq.gz")
    params:
        extra = "dedupe optical -Xmx4g -da" # suppress assertions
    resources:
        runtime = lambda wildcards, attempt: 90 + (attempt * 30),
        mem_mb = 16000
    log: 
        "output/{run}/log/clumpify.log"
    wrapper:
        WRAPPER_PREFIX + "master/bbtools/clumpify"


rule filterbytile:
    input:
        input = rules.clumpify.output.out
    output:
        out = temp("output/{run}/filterbytile.fq.gz")
    params:
        extra = "-Xmx4g -da"
    resources:
        runtime = 120,
        mem_mb = 16000
    log: 
        "output/{run}/log/filterbytile.log"
    wrapper:
        WRAPPER_PREFIX + "master/bbtools/filterbytile"


rule trim:
    input:
        input = rules.filterbytile.output.out
    output:
        out = temp("output/{run}/trimmed.fq.gz")
    params:
        extra = "ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=70 ref=adapters ftm=5 ordered -Xmx4g -da"
    resources:
        runtime = lambda wildcards, attempt: 90 + (attempt * 30),
        mem_mb = 16000
    log: 
        "output/{run}/log/trim.log"
    wrapper:
        WRAPPER_PREFIX + "master/bbtools/bbduk"


rule artifacts:
    input:
        input = rules.trim.output.out
    output:
        out = "output/{run}/filtered.fq.gz"
    params:
        extra = "k=31 ref=artifacts,phix ordered cardinality -Xmx4g -da"
    resources:
        runtime = lambda wildcards, attempt: 90 + (attempt * 20),
        mem_mb = 16000
    log: 
        "output/{run}/log/artifacts.log"
    wrapper:
        WRAPPER_PREFIX + "master/bbtools/bbduk"


# Remove host sequences
rule maphost:
    input:
        input = rules.artifacts.output.out,
        ref = HOST_GENOME
    output:
        outu = "output/{run}/unmaphost.fq.gz",
        outm = "output/{run}/maphost.fq.gz",
        statsfile = "output/{run}/maphost.txt"
    params:
        extra = "nodisk -Xmx24g"
    log: 
        "output/{run}/log/maphost.log"
    resources:
        runtime = lambda wildcards, attempt: 120 + (attempt * 60),
        mem_mb = 24000
    threads: 4
    wrapper:
        WRAPPER_PREFIX + "master/bbtools/bbwrap"


rule correct1:
    input:
        input = rules.maphost.output.outu
    output:
        out = temp("output/{run}/ecco.fq.gz")
    params:
        extra = "ecco mix vstrict ordered -Xmx4g -da"
    log: 
        "output/{run}/log/correct1.log"
    resources:
        runtime = lambda wildcards, attempt: 90 + (attempt * 30),
        mem_mb = 16000
    threads: 4
    wrapper:
        WRAPPER_PREFIX + "master/bbtools/bbmerge"


rule correct2:
    input:
        input = rules.correct1.output.out
    output:
        out = temp("output/{run}/eccc.fq.gz")
    params:
        extra = "passes=4 reorder -Xmx16g -da"
    log: 
        "output/{run}/log/correct2.log"
    resources:
        runtime = lambda wildcards, attempt: 90 + (attempt * 30),
        mem_mb = 16000
    wrapper:
        WRAPPER_PREFIX + "master/bbtools/clumpify"


rule correct3:
    input:
        input = rules.correct2.output.out
    output:
        out = temp("output/{run}/ecct.fq.gz")
    params:
        extra = "ecc k=62 ordered -Xmx16g -da"
    log: 
        "output/{run}/log/correct3.log"
    resources:
        runtime = lambda wildcards, attempt: 90 + (attempt * 30),
        mem_mb = 16000
    wrapper:
        WRAPPER_PREFIX + "master/bbtools/tadpole"


rule merge:
    input:
        input = rules.correct3.output.out
    output:
        out = temp("output/{run}/merged.fq.gz"),
        outu = temp("output/{run}/unmerged.fq.gz"),
        ihist = "output/{run}/ihist.txt"
    params:
        extra = "strict k=93 extend2=80 rem ordered -Xmx16g"
    log: 
        "output/{run}/log/merge.log"
    resources:
        runtime = lambda wildcards, attempt: 90 + (attempt * 30),
        mem_mb = 16000
    threads: 4
    wrapper:
        WRAPPER_PREFIX + "master/bbtools/bbmerge"


rule qtrim:
    input:
        input = rules.merge.output.outu
    output:
        out = temp("output/{run}/qtrimmed.fq.gz")
    params:
        extra = "qtrim=r trimq=10 minlen=70 ordered -Xmx4g"
    resources:
        runtime = lambda wildcards, attempt: 90 + (attempt * 30),
        mem_mb = 16000
    log: 
        "output/{run}/log/qtrim.log"
    wrapper:
        WRAPPER_PREFIX + "master/bbtools/bbduk"


rule concatenate:
    input:
      rules.merge.output.out, rules.qtrim.output.out
    output:
      temp("output/{run}/concatenated.fq.gz")
    resources:
        runtime = 120,
        mem_mb = 4000
    shell:
      "cat {input} > {output}"


rule assemble:
    input: 
      se = rules.concatenate.output[0]
    output: 
      contigs = "output/{run}/final.contigs.fa"
    params:
      extra = "--min-contig-len 1000"
    threads: 4
    log: 
        "output/{run}/log/assemble.log"
    shadow: 
      "minimal"
    resources:
        runtime = 2400,
        mem_mb = 96000
    wrapper:
      WRAPPER_PREFIX + "master/assembly/megahit"


# Calculate assembly coverage stats
# nodisk keeps index in memory, otherwise index will be written once to project root (ref/1) from first run to be processed 
# and reused for other unrelated runs.
# Key "input" will be parsed to "in", "input1" to "in1" etc.
rule coverage:
    input:
      ref = rules.assemble.output.contigs, 
      input = rules.concatenate.output[0]
    output:
      out = "output/{run}/final.contigs_aln.sam",
      covstats = "output/{run}/coverage.txt",
      statsfile = "output/{run}/mapcontigs.txt"
    params: 
      extra = "kfilter=22 subfilter=15 maxindel=80 nodisk"
    resources:
        runtime = lambda wildcards, attempt: 120 + (attempt * 60),
        mem_mb = 24000
    wrapper:
      WRAPPER_PREFIX + "master/bbtools/bbwrap"


# Tantan mask of low complexity DNA sequences
rule tantan:
    input:
        rules.assemble.output.contigs
    output:
        temp("output/{run}/tantan.fa")
    params:
        extra = "-x N" # mask low complexity using N
    resources:
        runtime = 120,
        mem_mb = 8000
    wrapper:
        WRAPPER_PREFIX + "master/tantan"


# Filter tantan output
# 1) Sequences > 50 nt of consecutive sequence without N
# 2) Sequences with >= 40% of total length of being masked
rule tantan_good:
    input:
        masked = rules.tantan.output[0]
    output:
        masked_filt = temp("output/{run}/tantangood.fa")
    params:
        min_length = 50,
        por_n = 40
    resources:
        runtime = 120
    wrapper:
        WRAPPER_PREFIX + "master/filter/masked"


# Split reads to smaller chunks for RepeatMasker
rule split_fasta:
    input:
        rules.tantan_good.output.masked_filt
    output:
        temp(expand("output/{{run}}/repeatmasker_{n}.fa", n = N))
    params:
        config["split_fasta"]["n_files"]
    resources:
        runtime = lambda wildcards, attempt: 90 + (attempt * 30),
        mem_mb = 4000
    wrapper:
        WRAPPER_PREFIX + "master/split-fasta"


# Repeatmasker
# Outputs are generated from input file names by RepeatMasker
# must have file extension '.masked'
# If no repetitive sequences were detected symlink output to input file
rule repeatmasker:
    input:
        fa = "output/{run}/repeatmasker_{n}.fa"
    output:
        masked = temp("output/{run}/repeatmasker_{n}.fa.masked"),
        out = temp("output/{run}/repeatmasker_{n}.fa.out"),
        cat = temp("output/{run}/repeatmasker_{n}.fa.cat"),
        tbl = "output/{run}/repeatmasker_{n}.fa.tbl"
    params:
        extra = "-qq"
    threads: 8
    resources:
        runtime = 1440,
        mem_mb = 16000
    singularity:
        "shub://tpall/repeatmasker-singularity"
    script:
        WRAPPER_PREFIX + "master/repeatmasker/wrapper.py"


# Filter repeatmasker output
# 1) Sequences > min_length nt of consecutive sequence without N
# 2) Sequences with >= % of total length of being masked
# input, output, and params names must match function arguments
rule repeatmasker_good:
    input:
        masked = rules.repeatmasker.output.masked,
        original = rules.repeatmasker.input.fa
    output:
        masked_filt = temp("output/{run}/repmaskedgood_{n}.fa"),
        original_filt = temp("output/{run}/unmaskedgood_{n}.fa")
    params:
        min_length = 100,
        por_n = 30
    resources:
        runtime = 120
    wrapper:
        WRAPPER_PREFIX + "master/filter/masked"

# Read QC stats
rule fastqc:
    input:
        rules.interleave.output.out
    output:
        html = "output/{run}/fastqc.html",
        zip = "output/{run}/fastqc.zip"
    resources:
        runtime = 120,
        mem_mb = 4000    
    wrapper:
        "0.27.1/bio/fastqc"


rule multiqc:
    input:
        "output/{run}/fastqc.zip",
        "output/{run}/maphost.txt",
        "output/{run}/coverage.txt",
        "output/{run}/mapcontigs.txt"
    output:
        report("output/{run}/multiqc.html", caption = "report/multiqc.rst", category = "Quality control")
    log:
        "output/{run}/log/multiqc.log"
    resources:
        runtime = 120,
        mem_mb = 4000    
    wrapper:
      WRAPPER_PREFIX + "master/multiqc"