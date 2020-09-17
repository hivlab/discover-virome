
def get_fastq(wildcards):
    fq_cols = [col for col in df.columns if "fq" in col]
    fqs = df.loc[(wildcards.group, wildcards.run), fq_cols].dropna()
    assert len(fq_cols) in [1, 2], "Enter one or two FASTQ file paths"
    if len(fq_cols) == 2:
        return {"in1": fqs[0], "in2": fqs[1]}
    else:
        return {"input": fqs[0]}
