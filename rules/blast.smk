
# Helper function to import tables
def safely_read_csv(path, **kwargs):
    try:
        return pd.read_csv(path, **kwargs)
    except pd.errors.EmptyDataError:
        pass


RANKS_OF_INTEREST = ["superkingdom", "order", "family", "genus", "species"]


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
      viruses = "output/blast/viruses.taxids",
      negative = "output/blast/negative.taxids"
    params: 
      viruses = 10239, 
      negative = [HOST_TAXID, 2, 12908]
    wrapper:
        wrapper_prefix + "master/blast/taxidslist"


# Blastn, megablast and blastx input, output, and params keys must match commandline blast option names. Please see https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a for all available options.
# Blast against nt virus database.
rule blastn_virus:
    input:
      query = "output/RM/{run}_repmaskedgood_{n}.fa",
      taxidlist = "output/blast/viruses.taxids"
    output:
      out = temp("output/blast/{run}_blastn-virus_{n}.tsv")
    params:
      program = "blastn",
      db = "nt_v5",
      evalue = 1e-4,
      max_hsps = 50,
      outfmt = "'6 qseqid sacc staxid pident length evalue'"
    threads: 8
    wrapper:
      BLAST_QUERY


# Filter blastn hits for the cutoff value.
rule parse_blastn_virus:
    input:
      query = "output/RM/{run}_repmaskedgood_{n}.fa",
      blast_result = rules.blastn_virus.output.out
    output:
      mapped = temp("output/blast/{run}_blastn-virus_{n}_mapped.tsv"),
      unmapped = temp("output/blast/{run}_blastn-virus_{n}_unmapped.fa")
    params:
      e_cutoff = 1e-5,
      outfmt = rules.blastn_virus.params.outfmt
    wrapper:
      PARSE_BLAST


# Blastx unmapped reads against nr virus database.
rule blastx_virus:
    input:
      query = rules.parse_blastn_virus.output.unmapped,
      taxidlist = "blast/viruses.taxids"
    output:
      out = temp("output/blast/{run}_blastx-virus_{n}.tsv")
    params:
      program = "blastx",
      task = "Blastx-fast",
      db = "nr_v5",
      evalue = 1e-2,
      db_soft_mask = 100,
      max_hsps = 50,
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
      mapped = temp("output/blast/{run}_blastx-virus_{n}_mapped.tsv"),
      unmapped = temp("output/blast/{run}_blastx-virus_{n}_unmapped.fa")
    params:
      e_cutoff = 1e-3,
      outfmt = rules.blastn_virus.params.outfmt
    wrapper:
      PARSE_BLAST


# Megablast against nt database.
rule megablast_nt:
    input:
      query = rules.parse_blastx_virus.output.unmapped if config["run_blastx"] else rules.parse_blastn_virus.output.unmapped,
      negative_taxidlist = "output/blast/negative.taxids"
    output:
      out = temp("output/blast/{run}_megablast-nt_{n}.tsv")
    params:
      program = "blastn",
      db = "nt_v5",
      task = "megablast",
      evalue = 1e-8,
      word_size = 16,
      max_hsps = 50,
      outfmt = rules.blastn_virus.params.outfmt
    threads: 8
    wrapper:
      BLAST_QUERY


# Filter megablast hits for the cutoff value.
rule parse_megablast_nt:
    input:
      query = rules.megablast_nt.input.query,
      blast_result = rules.megablast_nt.output.out
    output:
      mapped = temp("output/blast/{run}_megablast-nt_{n}_mapped.tsv"),
      unmapped = temp("output/blast/{run}_megablast-nt_{n}_unmapped.fa")
    params:
      e_cutoff = 1e-10,
      outfmt = rules.blastn_virus.params.outfmt
    wrapper:
      PARSE_BLAST


# Blastn against nt database.
rule blastn_nt:
    input:
      query = rules.parse_megablast_nt.output.unmapped,
      negative_taxidlist = "output/blast/negative.taxids"
    output:
      out = temp("output/blast/{run}_blastn-nt_{n}.tsv")
    params:
      program = "blastn",
      db = "nt_v5",
      task = "blastn",
      evalue = 1e-8,
      max_hsps = 50,
      outfmt = rules.blastn_virus.params.outfmt
    threads: 8
    wrapper:
      BLAST_QUERY


# Filter blastn records for the cutoff value.
rule parse_blastn_nt:
    input:
      query = rules.blastn_nt.input.query,
      blast_result = rules.blastn_nt.output.out
    output:
      mapped = temp("output/blast/{run}_blastn-nt_{n}_mapped.tsv"),
      unmapped = temp("output/blast/{run}_blastn-nt_{n}_unmapped.fa")
    params:
      e_cutoff = 1e-10,
      outfmt = rules.blastn_virus.params.outfmt
    wrapper:
      PARSE_BLAST


# Blastx unmapped sequences against nr database.
rule blastx_nr:
    input:
      query = rules.parse_blastn_nt.output.unmapped,
      negative_taxidlist = "blast/negative.taxids"
    output:
      out = temp("output/blast/{run}_blastx-nr_{n}.tsv")
    params:
      program = "blastx",
      task = "Blastx-fast",
      db = "nr_v5",
      evalue = 1e-2,
      max_hsps = 50,
      outfmt = rules.blastn_virus.params.outfmt
    threads: 8
    wrapper:
      BLAST_QUERY


# Filter blastx records for the cutoff value.
rule parse_blastx_nr:
    input:
      query = rules.blastx_nr.input.query,
      blast_result = rules.blastx_nr.output.out
    output:
      mapped = temp("output/blast/{run}_blastx-nr_{n}_mapped.tsv"),
      unmapped = temp("output/blast/{run}_blastx-nr_{n}_unmapped.fa")
    params:
      e_cutoff = 1e-3,
      outfmt = rules.blastn_virus.params.outfmt
    wrapper:
      PARSE_BLAST


# Filter sequences by division id.
# Saves hits with division id
rule classify_all:
  input:
    expand("output/blast/{{run}}_{blastresult}_{{n}}_mapped.tsv", blastresult = BLASTNR)
  output:
    temp("output/results/{run}_all_{n}.csv")
  params:
    pp_sway = 1, 
    ranks_of_interest = RANKS_OF_INTEREST,
    dbfile = TAXON_DB
  wrapper:
    BLAST_TAXONOMY


# Split classification rule outputs into viruses and non-viral
rule filter_viruses:
  input:
    expand("output/results/{{run}}_all_{n}.csv", n = N)
  output:
    viral = "output/results/{run}_viruses.csv",
    non_viral = "output/results/{run}_non-viral.csv"
  params:
    ranks = RANKS_OF_INTEREST
  run:
    tab = concatenate_tables(input, sep = ",", cols_to_integer = params.ranks)
    vir = tab[tab.superkingdom == 10239]
    non_vir = tab[tab.superkingdom != 10239]
    vir.to_csv(output.viral, index = False)
    non_vir.to_csv(output.non_viral, index = False)


# Merge unassigned sequences
rule merge_unassigned:
  input:
    expand("output/blast/{{run}}_blast{type}_{n}_unmapped.fa", type = "x-nr" if config["run_blastx"] else "n-nt", n = N)
  output:
    "output/results/{run}_unassigned.fa"
  shell:
    "cat {input} > {output}"


# Collect stats.
rule blast_stats:
  input:
    expand("output/blast/{{run}}_{blastresult}_{n}_unmapped.fa", blastresult = BLAST, n = N)
  output:
    "output/stats/{run}_blast.tsv"
  params:
    extra = "-T"
  wrapper:
    SEQ_STATS
