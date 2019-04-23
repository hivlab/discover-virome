
# Prepare taxonomy annotation tables.
rule prepare_taxonomy_data:
  input: config["names"], config["nodes"], config["division"]
  output:
      expand("taxonomy/{file}.csv", file = ["names", "nodes", "division"])
  conda:
    "../envs/tidyverse.yaml"
  script:
    "../scripts/prepare_taxonomy_data.R"

# Blastn, megablast and blastx input, output, and params keys must match commandline blast option names. Please see https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a for all available options.
# Blast against nt virus database.
rule blastn_virus:
    input:
      query = rules.repeatmasker_good.output.masked_filt
    output:
      out = "assemble/blast/{run}_blastn_virus.tsv"
    params:
      db = config["virus_nt"],
      task = "blastn",
      evalue = config["blastn_virus"]["evalue"],
      db_soft_mask = config["blastn_virus"]["db_soft_mask"],
      max_hsps = config["blastn_virus"]["max_hsps"],
      show_gis = True,
      num_threads = 8,
      outfmt = "'6 qseqid sgi pident length mismatch gapopen qstart qend sstart send evalue bitscore'"
    wrapper:
      config["wrappers"]["blast"]

# Filter blastn hits for the cutoff value.
rule parse_blastn_virus:
    input:
      query = rules.repeatmasker_good.output.masked_filt,
      blast_result = rules.blastn_virus.output.out
    output:
      mapped = "assemble/blast/{run}_blastn_virus_known-viral.tsv",
      unmapped = "assemble/blast/{run}_blastn_virus_unmapped.fa"
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
      out = "assemble/blast/{run}_blastx_virus.tsv"
    params:
      db = config["virus_nr"],
      word_size = 6,
      evalue = config["blastx_virus"]["evalue"],
      db_soft_mask = config["blastx_virus"]["db_soft_mask"],
      max_hsps = config["blastx_virus"]["max_hsps"],
      show_gis = True,
      num_threads = 8,
      outfmt = rules.blastn_virus.params.outfmt
    wrapper:
      config["wrappers"]["blast"]

# Filter blastn hits for the cutoff value.
rule parse_blastx_virus:
    input:
      query = rules.blastx_virus.input.query,
      blast_result = rules.blastx_virus.output.out
    output:
      mapped = "assemble/blast/{run}_blastx_virus_known-viral.tsv",
      unmapped = "assemble/blast/{run}_blastx_virus_unmapped.fa"
    params:
      e_cutoff = 1e-3,
      outfmt = rules.blastn_virus.params.outfmt
    wrapper:
      config["wrappers"]["parse_blast"]

# Filter sequences by division id.
# Saves hits with division id
rule classify_phages:
  input:
    [rules.parse_blastn_virus.output.mapped,
    rules.parse_blastx_virus.output.mapped] if config["run_blastx"] else rules.parse_blastn_virus.output.mapped,
    nodes = "taxonomy/nodes.csv"
  output:
    division = "assemble/results/{run}_phages.csv",
    other = "assemble/blast/{run}_candidate_viruses.csv"
  params:
    taxdb = config["vhunter"],
    division_id = 3
  wrapper:
    config["wrappers"]["blast_taxonomy"]

# Filter unmasked candidate virus reads.
rule unmasked_other:
    input:
      rules.classify_phages.output.other,
      rules.repeatmasker_good.output.original_filt
    output:
      "assemble/blast/{run}_candidate_viruses_unmasked.fa"
    conda:
      "../envs/biopython.yaml"
    script:
      "../scripts/unmasked_viral.py"

# Map reads against bacterial genomes.
rule bwa_mem_refbac:
  input:
    reads = [rules.unmasked_other.output]
  output:
    temp("assemble/blast/{run}_bac_mapped.sam")
  params:
    index = config["ref_bacteria"],
    sort = "samtools",
    sort_order = "queryname"
  log:
    "logs/{run}_bwa_map_refbac.log"
  threads: 2
  wrapper:
    "0.32.0/bio/bwa/mem"

# Extract unmapped reads.
rule unmapped_refbak:
  input:
    rules.bwa_mem_refbac.output
  output:
    "assemble/blast/{run}_refgenome_unmapped.bam"
  params:
    "-b -f 4" # bam output
  wrapper:
    "0.32.0/bio/samtools/view"

# Convert bam file to fastq file.
rule unmapped_refbaktofastq:
  input:
    rules.unmapped_refbak.output
  output:
    temp("assemble/blast/{run}_bac_unmapped.fq")
  log: 
    "logs/{run}_bamtofastq.log"
  wrapper:
    "https://bitbucket.org/tpall/snakemake-wrappers/raw/8e23fd260cdbed02450a7eb1796dce984d2e1f8f/bio/bedtools/bamtofastq"

# Convert fastq file to fasta file.
rule unmapped_refbaktofasta:
  input:
    rules.unmapped_refbaktofastq.output
  output:
    "assemble/blast/{run}_bac_unmapped.fa"
  shell:
    "cat {input} | sed -n '1~4s/^@/>/p;2~4p' > {output}"

# Subset repeatmasker masked reads using unmapped reads.
rule refbac_unmapped_masked:
    input:
      rules.unmapped_refbaktofasta.output,
      rules.repeatmasker_good.output.masked_filt
    output:
      temp("assemble/blast/{run}_bac_unmapped_masked.fa")
    conda:
      "../envs/biopython.yaml"
    script:
      "../scripts/unmapped_masked_ids.py"

# Megablast against nt database.
rule megablast_nt:
    input:
      query = rules.refbac_unmapped_masked.output
    output:
      out = "assemble/blast/{run}_megablast_nt.tsv"
    params:
      db = config["nt"],
      task = "megablast",
      evalue = config["megablast_nt"]["evalue"],
      word_size = config["megablast_nt"]["word_size"],
      max_hsps = config["megablast_nt"]["max_hsps"],
      show_gis = True,
      num_threads = 8,
      outfmt = rules.blastn_virus.params.outfmt
    wrapper:
      config["wrappers"]["blast"]

# Filter megablast hits for the cutoff value.
rule parse_megablast_nt:
    input:
      query = rules.refbac_unmapped_masked.output,
      blast_result = rules.megablast_nt.output.out
    output:
      mapped = "assemble/blast/{run}_megablast_nt_mapped.tsv",
      unmapped = "assemble/blast/{run}_megablast_nt_unmapped.fa"
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
      out = "assemble/blast/{run}_blastn_nt.tsv"
    params:
      db = config["nt"],
      task = "blastn",
      evalue = config["blastn_nt"]["evalue"],
      max_hsps = config["blastn_nt"]["max_hsps"],
      show_gis = True,
      num_threads = 8,
      outfmt = rules.blastn_virus.params.outfmt
    wrapper:
      config["wrappers"]["blast"]

# Filter blastn records for the cutoff value.
rule parse_blastn_nt:
    input:
      query = rules.blastn_nt.input.query,
      blast_result = rules.blastn_nt.output.out
    output:
      mapped = "assemble/blast/{run}_blastn_nt_mapped.tsv",
      unmapped = "assemble/blast/{run}_blastn_nt_unmapped.fa" if config["run_blastx"] else "assemble/results/{run}_unassigned.fa"
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
      out = "assemble/blast/{run}_blastx_nr.tsv"
    params:
      db = config["nr"],
      word_size = 6,
      evalue = config["blastx_nr"]["evalue"],
      max_hsps = config["blastx_nr"]["max_hsps"],
      show_gis = True,
      num_threads = 8,
      outfmt = rules.blastn_virus.params.outfmt
    wrapper:
      config["wrappers"]["blast"]

# Filter blastx records for the cutoff value.
rule parse_blastx_nr:
    input:
      query = rules.blastx_nr.input.query,
      blast_result = rules.blastx_nr.output.out
    output:
      mapped = "assemble/blast/{run}_blastx_nr_mapped.tsv",
      unmapped = "assemble/results/{run}_unassigned.fa" if config["run_blastx"] else "{run}_None"
    params:
      e_cutoff = 1e-3,
      outfmt = rules.blastn_virus.params.outfmt
    wrapper:
      config["wrappers"]["parse_blast"]

# Filter sequences by division id.
# Saves hits with division id
rule classify_phages_viruses:
  input:
    [rules.parse_megablast_nt.output.mapped, rules.parse_blastn_nt.output.mapped, rules.parse_blastx_nr.output.mapped] if config["run_blastx"] else [rules.parse_megablast_nt.output.mapped, rules.parse_blastn_nt.output.mapped],
    nodes = "taxonomy/nodes.csv"
  output:
    division = "assemble/results/{run}_phages_viruses.csv",
    other = "assemble/results/{run}_non_viral.csv"
  params:
    taxdb = config["vhunter"],
    division_id = [3, 9] # pool phages and viruses
  wrapper:
    config["wrappers"]["blast_taxonomy"]

if config["zenodo"]["deposition_id"]:
  rule upload:
    input: 
      "assemble/results/{run}_{result}.{ext}"
    output: 
      ZEN.remote(expand("{deposition_id}/files/assemble/results/{{run, [^_]+}}_{{result}}.{{ext}}", deposition_id = config["zenodo"]["deposition_id"]))
    shell: 
      "cp {input} {output}"
