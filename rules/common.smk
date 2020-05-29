
# Check file exists
def file_exists(file):
    try:
        with open(file, 'r') as fh:
            print("{} is set up correctly".format(file))
    except FileNotFoundError:
        ("Could not find {}").format(file)

# Helper functions to get input files
def get_fastq(wildcards):
    """Get fraction read file paths from samples.tsv"""
    urls = RUNS.loc[wildcards.run, ['fq1', 'fq2']]
    return list(urls)

def get_frac(wildcards):
    """Get fraction of reads to be sampled from samples.tsv"""
    frac = RUNS.loc[wildcards.run, ['frac']][0]
    return frac

# Helper function to import tables
def safely_read_csv(path, **kwargs):
    try:
        return pd.read_csv(path, **kwargs)
    except pd.errors.EmptyDataError:
        pass

# Helper function to concatenate output tables
def concatenate_tables(input, sep="\s+", cols_to_integer=None):
    frames = [safely_read_csv(f, sep=sep) for f in input]
    frames_concatenated = pd.concat(frames, keys=input, sort=False)
    if cols_to_integer:
        frames_concatenated[cols_to_integer] = frames_concatenated[
            cols_to_integer
        ].apply(lambda x: pd.Series(x, dtype="Int64"))
    return frames_concatenated
