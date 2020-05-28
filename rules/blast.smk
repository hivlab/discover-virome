
# Helper function to import tables
def safely_read_csv(path, **kwargs):
    try:
        return pd.read_csv(path, **kwargs)
    except pd.errors.EmptyDataError:
        pass


RANKS_OF_INTEREST = ["superkingdom", "order", "family", "genus", "species"]
VIRUSES_TAXID = 10239

def concatenate_tables(input, sep="\s+", cols_to_integer=None):
    frames = [safely_read_csv(f, sep=sep) for f in input]
    frames_concatenated = pd.concat(frames, keys=input, sort=False)
    if cols_to_integer:
        frames_concatenated[cols_to_integer] = frames_concatenated[
            cols_to_integer
        ].apply(lambda x: pd.Series(x, dtype="Int64"))
    return frames_concatenated


# Creates to required outputs viruses.taxids and negative.taxids.
# Output directory can be changed.
# Additional negative taxids (all listed taxids except viruses) can be added via params.
# Shadow=full ensures that only required outputs will be saved. 
rule taxids_list:
    output:
      viruses = "output/blast/viruses.taxids"
    params: 
      viruses = VIRUSES_TAXID
    wrapper:
        WRAPPER_PREFIX + "master/blast/taxidslist"


# Blast input, output, and params keys must match commandline blast option names. 
# Please see https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a 
# for all available options.
# Blast against nt virus database.
rule blastn_virus:
    input:
      query = "output/{run}/repmaskedgood_{n}.fa",
      taxidlist = "output/blast/viruses.taxids"
    output:
      out = temp("output/{run}/blastn-virus_{n}.tsv")
    params:
      program = "blastn",
      db = "nt_v5",
      evalue = 1e-6,
      max_hsps = 1,
      outfmt = "'6 qseqid sacc staxid pident length evalue'"
    threads: 8
    wrapper:
      BLAST_QUERY


# Filter blastn hits for the cutoff value.
rule parse_blastn_virus:
    input:
      query = "output/{run}/repmaskedgood_{n}.fa",
      blast_result = rules.blastn_virus.output.out
    output:
      mapped = temp("output/{run}/blastn-virus_{n}_mapped.tsv"),
      unmapped = temp("output/{run}/blastn-virus_{n}_unmapped.fa")
    params:
      e_cutoff = 1e-6,
      outfmt = rules.blastn_virus.params.outfmt
    wrapper:
      PARSE_BLAST


# Blastx unmapped reads against nr virus database.
rule blastx_virus:
    input:
      query = rules.parse_blastn_virus.output.unmapped,
      taxidlist = "output/blast/viruses.taxids"
    output:
      out = temp("output/{run}/blastx-virus_{n}.tsv")
    params:
      program = "blastx",
      task = "blastx-fast",
      db = "nr_v5",
      evalue = 1e-2,
      max_hsps = 1,
      outfmt = rules.blastn_virus.params.outfmt
    threads: 8
    wrapper:
      BLAST_QUERY


# Filter blastn hits for the cutoff value.
rule parse_blastx_virus:
    input:
      query = rules.blastx_virus.input.query,
      blast_result = rules.blastx_virus.output.out
    output:
      mapped = temp("output/{run}/blastx-virus_{n}_mapped.tsv"),
      unmapped = temp("output/{run}/blastx-virus_{n}_unmapped.fa")
    params:
      e_cutoff = 1e-3,
      outfmt = rules.blastn_virus.params.outfmt
    wrapper:
      PARSE_BLAST


# Filter sequences by division id.
# Saves hits with division id
rule classify_all:
  input:
    expand("output/{{run}}/{blastresult}_{{n}}_mapped.tsv", blastresult = BLASTV)
  output:
    temp("output/{run}/all_{n}.csv")
  params:
    pp_sway = 1, 
    ranks_of_interest = RANKS_OF_INTEREST,
    dbfile = TAXON_DB
  wrapper:
    BLAST_TAXONOMY


# Split classification rule outputs into viruses and non-viral
rule filter_viruses:
  input:
    expand("output/{{run}}/all_{n}.csv", n = N)
  output:
    viral = "output/{run}/viruses.csv",
    non_viral = "output/{run}/non-viral.csv"
  params:
    ranks = RANKS_OF_INTEREST
  run:
    tab = concatenate_tables(input, sep = ",", cols_to_integer = params.ranks)
    vir = tab[tab.superkingdom == VIRUSES_TAXID]
    non_vir = tab[tab.superkingdom != VIRUSES_TAXID]
    vir.to_csv(output.viral, index = False)
    non_vir.to_csv(output.non_viral, index = False)


# Merge unassigned sequences
rule merge_unassigned:
  input:
    expand("output/{{run}}/blast{type}_{n}_unmapped.fa", type = "x-virus" if config["run_blastx"] else "n-virus", n = N)
  output:
    "output/{run}/unassigned.fa"
  shell:
    "cat {input} > {output}"

