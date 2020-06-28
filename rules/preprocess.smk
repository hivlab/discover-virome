
FTP = FTPRemoteProvider(username = config["username"], password = config["password"])

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
        extra = lambda wildcards, resources: f"-Xmx{resources.mem_mb / 1000:.0f}g"
    resources:
        runtime = 120,
        mem_mb = 4000
    log:
        "output/{run}/log/interleave.txt"
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/reformat"


# Remove PCR and optical duplicates
rule clumpify:
    input:
        input = rules.interleave.output.out
    output:
        out = temp("output/{run}/clumpify.fq.gz")
    params:
        extra = lambda wildcards, resources: f"dedupe optical -Xmx{resources.mem_mb / 1000:.0f}g -da"
    resources:
        runtime = lambda wildcards, attempt: 90 + (attempt * 30),
        mem_mb = 16000
    log: 
        "output/{run}/log/clumpify.log"
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/clumpify"


rule filterbytile:
    input:
        input = rules.clumpify.output.out
    output:
        out = temp("output/{run}/filterbytile.fq.gz")
    params:
        extra = lambda wildcards, resources: f"-Xmx{resources.mem_mb / 1000:.0f}g -da"
    resources:
        runtime = 120,
        mem_mb = 16000
    log: 
        "output/{run}/log/filterbytile.log"
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/filterbytile"


rule trim:
    input:
        input = rules.filterbytile.output.out
    output:
        out = temp("output/{run}/trimmed.fq.gz")
    params:
        extra = lambda wildcards, resources: f"ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=70 ref=adapters ftm=5 ordered -Xmx{resources.mem_mb / 1000:.0f}g -da"
    resources:
        runtime = lambda wildcards, attempt: 90 + (attempt * 30),
        mem_mb = 16000
    log: 
        "output/{run}/log/trim.log"
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/bbduk"


rule artifacts:
    input:
        input = rules.trim.output.out
    output:
        out = "output/{run}/filtered.fq.gz"
    params:
        extra = lambda wildcards, resources: f"k=31 ref=artifacts,phix ordered cardinality -Xmx{resources.mem_mb / 1000:.0f}g -da"
    resources:
        runtime = lambda wildcards, attempt: 90 + (attempt * 20),
        mem_mb = 16000
    log: 
        "output/{run}/log/artifacts.log"
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/bbduk"


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
        extra = lambda wildcards, resources: f"nodisk -Xmx{resources.mem_mb / 1000:.0f}g"
    log: 
        "output/{run}/log/maphost.log"
    resources:
        runtime = lambda wildcards, attempt: 90 + (attempt * 30),
        mem_mb = 80000
    threads: 8
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/bbwrap"


rule correct1:
    input:
        input = rules.maphost.output.outu
    output:
        out = temp("output/{run}/ecco.fq.gz")
    params:
        extra = lambda wildcards, resources: f"ecco mix vstrict ordered -Xmx{resources.mem_mb / 1000:.0f}g -da"
    log: 
        "output/{run}/log/correct1.log"
    resources:
        runtime = lambda wildcards, attempt: 90 + (attempt * 30),
        mem_mb = 8000
    threads: 8
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/bbmerge"


rule correct2:
    input:
        input = rules.correct1.output.out
    output:
        out = temp("output/{run}/eccc.fq.gz")
    params:
        extra = lambda wildcards, resources: f"passes=4 reorder -Xmx{resources.mem_mb / 1000:.0f}g -da"
    log: 
        "output/{run}/log/correct2.log"
    resources:
        runtime = lambda wildcards, attempt: 90 + (attempt * 30),
        mem_mb = 32000
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/clumpify"


rule correct3:
    input:
        input = rules.correct2.output.out
    output:
        out = temp("output/{run}/ecct.fq.gz")
    params:
        extra = lambda wildcards, resources: f"ecc k=62 ordered -Xmx{resources.mem_mb / 1000:.0f}g -da"
    log: 
        "output/{run}/log/correct3.log"
    resources:
        runtime = lambda wildcards, attempt: 90 + (attempt * 30),
        mem_mb = 32000
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/tadpole"


rule merge:
    input:
        input = rules.correct3.output.out
    output:
        out = temp("output/{run}/merged.fq.gz"),
        outu = temp("output/{run}/unmerged.fq.gz"),
        ihist = "output/{run}/ihist.txt"
    params:
        extra = lambda wildcards, resources: f"strict k=93 extend2=80 rem ordered -Xmx{resources.mem_mb / 1000:.0f}g"
    log: 
        "output/{run}/log/merge.log"
    resources:
        runtime = lambda wildcards, attempt: 90 + (attempt * 30),
        mem_mb = 32000
    threads: 8
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/bbmerge"


rule qtrim:
    input:
        input = rules.merge.output.outu
    output:
        out = temp("output/{run}/qtrimmed.fq.gz")
    params:
        extra = lambda wildcards, resources: f"qtrim=r trimq=10 minlen=70 ordered -Xmx{resources.mem_mb / 1000:.0f}g"
    resources:
        runtime = lambda wildcards, attempt: 90 + (attempt * 30),
        mem_mb = 16000
    log: 
        "output/{run}/log/qtrim.log"
    wrapper:
        f"{WRAPPER_PREFIX}/master/bbtools/bbduk"

