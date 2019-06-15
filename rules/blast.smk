
# Prepare taxonomy annotation tables.
rule prepare_taxonomy_data:
  input: config["names"], config["nodes"], config["division"]
  output:
      expand("taxonomy/{file}.csv", file = ["names", "nodes", "division"])
  wrapper:
    "https://raw.githubusercontent.com/avilab/vs-wrappers/master/prepare_taxonomy_data"

# Blastn, megablast and blastx input, output, and params keys must match commandline blast option names. Please see https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a for all available options.
# Blast against nt virus database.
rule blastn_virus:
    input:
      query = "assemble/mask/{run}_repmaskedgood_{n}.fa"
    output:
      out = temp("assemble/blast/{run}_blastn_virus_{n}.tsv")
    params:
      db = config["virus_nt"],
      task = "blastn",
      evalue = config["blastn_virus"]["evalue"],
      db_soft_mask = config["blastn_virus"]["db_soft_mask"],
      max_hsps = config["blastn_virus"]["max_hsps"],
      show_gis = True,
      num_threads = 2,
      outfmt = "'6 qseqid sgi pident length mismatch gapopen qstart qend sstart send evalue bitscore'"
    threads: 2
    wrapper:
      config["wrappers"]["blast"]

# Filter blastn hits for the cutoff value.
rule parse_blastn_virus:
    input:
      query = "assemble/mask/{run}_repmaskedgood_{n}.fa",
      blast_result = rules.blastn_virus.output.out
    output:
      mapped = temp("assemble/blast/{run}_blastn_virus_mapped_{n}.tsv"),
      unmapped = temp("assemble/blast/{run}_blastn_virus_unmapped_{n}.fa")
    params:
      e_cutoff = 1e-5,
      outfmt = rules.blastn_virus.params.outfmt
    wrapper:
      config["wrappers"]["parse_blast"]

# Blastx unmapped reads against nr virus database.
rule blastx_virus:
    input:
      query = rules.parse_blastn_virus.output.unmapped
    output:
      out = temp("assemble/blast/{run}_blastx_virus_{n}.tsv")
    params:
      db = config["virus_nr"],
      word_size = 6,
      evalue = config["blastx_virus"]["evalue"],
      db_soft_mask = config["blastx_virus"]["db_soft_mask"],
      max_hsps = config["blastx_virus"]["max_hsps"],
      show_gis = True,
      num_threads = 2,
      outfmt = rules.blastn_virus.params.outfmt
    threads: 2
    wrapper:
      config["wrappers"]["blast"]

# Filter blastn hits for the cutoff value.
rule parse_blastx_virus:
    input:
      query = rules.blastx_virus.input.query,
      blast_result = rules.blastx_virus.output.out
    output:
      mapped = temp("assemble/blast/{run}_blastx_virus_mapped_{n}.tsv"),
      unmapped = temp("assemble/blast/{run}_blastx_virus_unmapped_{n}.fa")
    params:
      e_cutoff = 1e-3,
      outfmt = rules.blastn_virus.params.outfmt
    wrapper:
      config["wrappers"]["parse_blast"]

if config["run_blastx"]:
  rule merge_virus_fasta:
    input:
      rules.parse_blastn_virus.output.unmapped, rules.parse_blastx_virus.output.unmapped
    output:
      temp("assemble/blast/{run}_blast_virus_unmapped_{n}.fa")
    shell:
      "cat {input} > {output}"

# Megablast against nt database.
rule megablast_nt:
    input:
      query = rules.merge_virus_fasta.output if config["run_blastx"] else rules.parse_blastn_virus.output.unmapped
    output:
      out = temp("assemble/blast/{run}_megablast_nt_{n}.tsv")
    params:
      db = config["nt"],
      task = "megablast",
      evalue = config["megablast_nt"]["evalue"],
      word_size = config["megablast_nt"]["word_size"],
      max_hsps = config["megablast_nt"]["max_hsps"],
      show_gis = True,
      num_threads = 2,
      outfmt = rules.blastn_virus.params.outfmt
    threads: 2
    wrapper:
      config["wrappers"]["blast"]

# Filter megablast hits for the cutoff value.
rule parse_megablast_nt:
    input:
      query = rules.megablast_nt.input.query,
      blast_result = rules.megablast_nt.output.out
    output:
      mapped = temp("assemble/blast/{run}_megablast_nt_mapped_{n}.tsv"),
      unmapped = temp("assemble/blast/{run}_megablast_nt_unmapped_{n}.fa")
    params:
      e_cutoff = 1e-10,
      outfmt = rules.blastn_virus.params.outfmt
    wrapper:
      config["wrappers"]["parse_blast"]

# Blastn against nt database.
rule blastn_nt:
    input:
      query = rules.parse_megablast_nt.output.unmapped
    output:
      out = temp("assemble/blast/{run}_blastn_nt_{n}.tsv")
    params:
      db = config["nt"],
      task = "blastn",
      evalue = config["blastn_nt"]["evalue"],
      max_hsps = config["blastn_nt"]["max_hsps"],
      show_gis = True,
      num_threads = 2,
      outfmt = rules.blastn_virus.params.outfmt
    threads: 2
    wrapper:
      config["wrappers"]["blast"]

# Filter blastn records for the cutoff value.
rule parse_blastn_nt:
    input:
      query = rules.blastn_nt.input.query,
      blast_result = rules.blastn_nt.output.out
    output:
      mapped = temp("assemble/blast/{run}_blastn_nt_mapped_{n}.tsv"),
      unmapped = temp("assemble/blast/{run}_blastn_nt_unmapped_{n}.fa") if config["run_blastx"] else temp("assemble/blast/{run}_unassigned_{n}.fa")
    params:
      e_cutoff = 1e-10,
      outfmt = rules.blastn_virus.params.outfmt
    wrapper:
      config["wrappers"]["parse_blast"]

# Blastx unmapped sequences against nr database.
rule blastx_nr:
    input:
      query = rules.parse_blastn_nt.output.unmapped
    output:
      out = temp("assemble/blast/{run}_blastx_nr_{n}.tsv")
    params:
      db = config["nr"],
      word_size = 6,
      evalue = config["blastx_nr"]["evalue"],
      max_hsps = config["blastx_nr"]["max_hsps"],
      show_gis = True,
      num_threads = 2,
      outfmt = rules.blastn_virus.params.outfmt
    threads: 2
    wrapper:
      config["wrappers"]["blast"]

# Filter blastx records for the cutoff value.
rule parse_blastx_nr:
    input:
      query = rules.blastx_nr.input.query,
      blast_result = rules.blastx_nr.output.out
    output:
      mapped = temp("assemble/blast/{run}_blastx_nr_mapped_{n}.tsv"),
      unmapped = temp("assemble/blast/{run}_unassigned_{n}.fa") if config["run_blastx"] else temp("{run}_None_{n}")
    params:
      e_cutoff = 1e-3,
      outfmt = rules.blastn_virus.params.outfmt
    wrapper:
      config["wrappers"]["parse_blast"]

# Merge blast results for classification
rule merge_blast_results:
  input: expand("assemble/blast/{{run}}_{{blastresult}}_mapped_{n}.tsv", n = N)
  output: "assemble/blast/{run}_{blastresult}_mapped.tsv"
  run:
    merged = pd.concat(input)
    with open(output, "w") as out:
      out.write(merged)

# Merge unassigned sequences
rule merge_unassigned:
  input: expand("assemble/blast/{{run}}_unassigned_{n}.fa", n = N)
  output: "assemble/results/{run}_unassigned.fa"
  shell:
    "cat {input} > {output}"

# Filter sequences by division id.
# Saves hits with division id
rule classify_phages_viruses:
  input:
    expand("assemble/blast/{{run}}_{blastresult}_mapped.tsv", blastresult = ["blastn_virus", "blastx_virus", "megablast_nt", "blastn_nt", "blastx_nr"] if config["run_blastx"] else ["blastn_virus", "megablast_nt", "blastn_nt"]),
    nodes = "taxonomy/nodes.csv"
  output:
    division = "assemble/results/{run}_phages_viruses.csv",
    other = "assemble/results/{run}_non_viral.csv"
  params:
    taxdb = config["vhunter"],
    division_id = [3, 9] # pool phages and viruses
  wrapper:
    config["wrappers"]["blast_taxonomy"]

# Assign unique taxons to blast queries
rule query_taxid:
  input:
    rules.classify_phages_viruses.output.division
  output:
    "assemble/results/{run}_query_taxid.csv"
  wrapper:
    "https://raw.githubusercontent.com/avilab/vs-wrappers/master/unique_taxons"

# Upload results to Zenodo.
if config["zenodo"]["deposition_id"]:
  rule upload_parsed:
    input:
      rules.query_taxid.output
    output:
      ZEN.remote(expand("{deposition_id}/files/assemble/results/{{run}}_query_taxid.csv", deposition_id = config["zenodo"]["deposition_id"]))
    group: "upload"
    shell:
      "cp {input} {output}"

  rule upload:
    input: 
      "assemble/results/{run}_{result}.{ext}"
    output: 
      ZEN.remote(expand("{deposition_id}/files/assemble/results/{{run, [^_]+}}_{{result}}.{{ext}}", deposition_id = config["zenodo"]["deposition_id"]))
    shell: 
      "cp {input} {output}"

# Merge unmapped seqs for stats
rule merge_blast_unmapped:
  input: expand("assemble/blast/{{run}}_{{blastresult}}_unmapped_{n}.fa", n = N)
  output: temp("assemble/blast/{run}_{blastresult}_unmapped.fa")
  shell:
    "cat {input} > {output}"

# Collect stats
rule blast_stats:
  input:
    ["assemble/blast/{run}_blastn_virus_unmapped.fa",
    "assemble/blast/{run}_blastx_virus_unmapped.fa",
    "assemble/blast/{run}_megablast_nt_unmapped.fa",
    "assemble/blast/{run}_blastn_nt_unmapped.fa",
    "assemble/results/{run}_unassigned.fa"] if config["run_blastx"] else ["assemble/blast/{run}_blastn_virus_unmapped.fa",
    "assemble/blast/{run}_megablast_nt_unmapped.fa",
    "assemble/results/{run}_unassigned.fa"]
  output:
    "assemble/stats/{run}_blast.tsv"
  params:
    extra = "-T"
  wrapper:
    config["wrappers"]["stats"]
