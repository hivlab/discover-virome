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
    log:
        "output/{group}/{run}/log/interleave.txt",
    params:
        extra="",
    resources:
        runtime=360,
        mem_mb=4000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.6/bbtools/reformat"


rule trim:
    input:
        input=rules.interleave.output.out,
    output:
        out=temp("output/{group}/{run}/trimmed.fq.gz"),
    log:
        "output/{group}/{run}/log/trim.log",
    params:
        extra=(
            "ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=70 ref=adapters ftm=5 ordered"
        ),
    resources:
        runtime=lambda wildcards, attempt: 180 + (attempt * 60),
        mem_mb=16000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.6/bbtools/bbduk"


rule artifacts:
    input:
        input=rules.trim.output.out,
    output:
        out=temp("output/{group}/{run}/filtered.fq.gz"),
    log:
        "output/{group}/{run}/log/artifacts.log",
    params:
        extra="k=31 ref=artifacts,phix ordered cardinality",
    resources:
        runtime=lambda wildcards, attempt: 180 + (attempt * 60),
        mem_mb=16000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.6/bbtools/bbduk"


# Remove host sequences
rule maphost:
    input:
        input=rules.artifacts.output.out,
        ref=HOST_GENOME,
    output:
        outu=temp("output/{group}/{run}/unmaphost.fq.gz"),
        outm=temp("output/{group}/{run}/maphost.fq.gz"),
        statsfile="output/{group}/{run}/maphost.txt",
    log:
        "output/{group}/{run}/log/maphost.log",
    params:
        extra="minratio=0.9 maxindel=3 bwr=0.16 bw=12 fast minhits=2 qtrim=r trimq=10 untrim idtag kfilter=25 maxsites=1 k=14 nodisk",
    resources:
        runtime=lambda wildcards, attempt: attempt * 360,
        mem_mb=80000,
    threads: 8
    wrapper:
        f"{WRAPPER_PREFIX}/v0.6/bbtools/bbwrap"


rule correct1:
    input:
        input=rules.maphost.output.outu,
    output:
        out=temp("output/{group}/{run}/ecco.fq.gz"),
    log:
        "output/{group}/{run}/log/correct1.log",
    params:
        extra="ecco mix vstrict ordered",
    resources:
        runtime=360,
        mem_mb=8000,
    threads: 8
    wrapper:
        f"{WRAPPER_PREFIX}/v0.6/bbtools/bbmerge"


rule correct2:
    input:
        input=rules.correct1.output.out,
    output:
        out=temp("output/{group}/{run}/ecct.fq.gz"),
    log:
        "output/{group}/{run}/log/correct3.log",
    params:
        extra="mode=correct k=50 ordered",
    resources:
        runtime=lambda wildcards, attempt: attempt * 720,
        mem_mb=lambda wildcards, input: round(32000 + 6 * input.size_mb),
    wrapper:
        f"{WRAPPER_PREFIX}/v0.6/bbtools/tadpole"


rule merge:
    input:
        input=rules.correct2.output.out,
    output:
        out=temp("output/{group}/{run}/merged.fq.gz"),
        outu=temp("output/{group}/{run}/unmerged.fq.gz"),
        ihist="output/{group}/{run}/ihist.txt",
    log:
        "output/{group}/{run}/log/merge.log",
    params:
        extra="strict k=93 extend2=80 rem ordered",
    resources:
        runtime=lambda wildcards, attempt: 180 + (attempt * 60),
        mem_mb=lambda wildcards, input: round(32000 + 6 * input.size_mb),
    threads: 8
    wrapper:
        f"{WRAPPER_PREFIX}/v0.6/bbtools/bbmerge"


rule qtrim:
    input:
        input=rules.merge.output.outu,
    output:
        out=temp("output/{group}/{run}/qtrimmed.fq.gz"),
    log:
        "output/{group}/{run}/log/qtrim.log",
    params:
        extra="maq=10 qtrim=r trimq=10 ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=100 ref=adapters ftm=5 ordered",
    resources:
        runtime=lambda wildcards, attempt: 180 + (attempt * 60),
        mem_mb=16000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/bbtools/bbduk"
