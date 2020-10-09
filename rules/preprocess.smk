# Convert reads to interleaved format
rule interleave:
    input:
        unpack(get_fastq),
    output:
        out=temp("output/{group}/{run}/interleaved.fq.gz"),
        bhist="output/{group}/{run}/bhist.txt",
        qhist="output/{group}/{run}/qhist.txt",
        aqhist="output/{group}/{run}/aqhist.txt",
        bqhist="output/{group}/{run}/bqhist.txt",
        lhist="output/{group}/{run}/lhist.txt",
        gchist="output/{group}/{run}/gchist.txt",
    params:
        extra=lambda wildcards, resources: f"-Xmx{resources.mem_mb / 1000:.0f}g",
    resources:
        runtime=120,
        mem_mb=4000,
    log:
        "output/{group}/{run}/log/interleave.txt",
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/bbtools/reformat"


# Remove PCR and optical duplicates
rule clumpify:
    input:
        input=rules.interleave.output.out,
    output:
        out=temp("output/{group}/{run}/clumpify.fq.gz"),
    params:
        extra=(
            lambda wildcards, resources: f"dedupe optical -Xmx{resources.mem_mb / 1000:.0f}g -da"
        ),
    resources:
        runtime=lambda wildcards, attempt: 90 + (attempt * 30),
        mem_mb=16000,
    log:
        "output/{group}/{run}/log/clumpify.log",
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/bbtools/clumpify"


rule filterbytile:
    input:
        input=rules.clumpify.output.out,
    output:
        out=temp("output/{group}/{run}/filterbytile.fq.gz"),
    params:
        extra=lambda wildcards, resources: f"-Xmx{resources.mem_mb / 1000:.0f}g -da",
    resources:
        runtime=120,
        mem_mb=16000,
    log:
        "output/{group}/{run}/log/filterbytile.log",
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/bbtools/filterbytile"


rule trim:
    input:
        input=rules.filterbytile.output.out,
    output:
        out=temp("output/{group}/{run}/trimmed.fq.gz"),
    params:
        extra=(
            lambda wildcards, resources: f"ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=70 ref=adapters ftm=5 ordered -Xmx{resources.mem_mb / 1000:.0f}g -da"
        ),
    resources:
        runtime=lambda wildcards, attempt: 90 + (attempt * 30),
        mem_mb=16000,
    log:
        "output/{group}/{run}/log/trim.log",
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/bbtools/bbduk"


rule artifacts:
    input:
        input=rules.trim.output.out,
    output:
        out=temp("output/{group}/{run}/filtered.fq.gz"),
    params:
        extra=(
            lambda wildcards, resources: f"k=31 ref=artifacts,phix ordered cardinality -Xmx{resources.mem_mb / 1000:.0f}g -da"
        ),
    resources:
        runtime=lambda wildcards, attempt: 90 + (attempt * 20),
        mem_mb=16000,
    log:
        "output/{group}/{run}/log/artifacts.log",
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/bbtools/bbduk"


# Remove host sequences
rule maphost:
    input:
        input=rules.artifacts.output.out,
        ref=HOST_GENOME,
    output:
        outu=temp("output/{group}/{run}/unmaphost.fq.gz"),
        outm=temp("output/{group}/{run}/maphost.fq.gz"),
        statsfile="output/{group}/{run}/maphost.txt",
    params:
        extra=lambda wildcards, resources: f"nodisk -Xmx{resources.mem_mb}m",
    log:
        "output/{group}/{run}/log/maphost.log",
    resources:
        runtime=lambda wildcards, attempt: attempt * 120,
        mem_mb=80000,
    threads: 8
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/bbtools/bbwrap"


rule correct1:
    input:
        input=rules.maphost.output.outu,
    output:
        out=temp("output/{group}/{run}/ecco.fq.gz"),
    params:
        extra=(
            lambda wildcards, resources: f"ecco mix vstrict ordered -Xmx{resources.mem_mb}m -da"
        ),
    log:
        "output/{group}/{run}/log/correct1.log",
    resources:
        runtime=120,
        mem_mb=8000,
    threads: 8
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/bbtools/bbmerge"


rule correct2:
    input:
        input=rules.correct1.output.out,
    output:
        out=temp("output/{group}/{run}/eccc.fq.gz"),
    params:
        extra=(
            lambda wildcards, resources: f"passes=4 reorder -Xmx{resources.mem_mb}m -da"
        ),
    log:
        "output/{group}/{run}/log/correct2.log",
    resources:
        runtime=lambda wildcards, attempt: attempt * 240,
        mem_mb=lambda wildcards, input: round(4000 + 3 * input.size_mb),
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/bbtools/clumpify"


rule correct3:
    input:
        input=rules.correct2.output.out,
    output:
        out=temp("output/{group}/{run}/ecct.fq.gz"),
    params:
        extra=(
            lambda wildcards, resources: f"mode=correct k=62 ordered -Xmx{resources.mem_mb}m -da"
        ),
    log:
        "output/{group}/{run}/log/correct3.log",
    resources:
        runtime=lambda wildcards, attempt: attempt * 120,
        mem_mb=lambda wildcards, input: round(32000 + 6 * input.size_mb),
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/bbtools/tadpole"


rule merge:
    input:
        input=rules.correct3.output.out,
    output:
        out=temp("output/{group}/{run}/merged.fq.gz"),
        outu=temp("output/{group}/{run}/unmerged.fq.gz"),
        ihist="output/{group}/{run}/ihist.txt",
    params:
        extra=(
            lambda wildcards, resources: f"strict k=93 extend2=80 rem ordered -Xmx{resources.mem_mb}m"
        ),
    log:
        "output/{group}/{run}/log/merge.log",
    resources:
        runtime=lambda wildcards, attempt: 90 + (attempt * 30),
        mem_mb=lambda wildcards, input: round(32000 + 6 * input.size_mb),
    threads: 8
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/bbtools/bbmerge"


rule qtrim:
    input:
        input=rules.merge.output.outu,
    output:
        out=temp("output/{group}/{run}/qtrimmed.fq.gz"),
    params:
        extra=(
            lambda wildcards, resources: f"maq=10 qtrim=r trimq=10 ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=100 ref=adapters ftm=5 ordered qin=33 ordered -Xmx{resources.mem_mb / 1000:.0f}g"
        ),
    resources:
        runtime=lambda wildcards, attempt: 90 + (attempt * 30),
        mem_mb=16000,
    log:
        "output/{group}/{run}/log/qtrim.log",
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/bbtools/bbduk"
