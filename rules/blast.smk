
RANKS_OF_INTEREST = ["superkingdom", "order", "family", "genus", "species"]
VIRUSES_TAXID = 10239


# Creates to required outputs viruses.taxids and negative.taxids.
# Output directory can be changed.
# Additional negative taxids (all listed taxids except viruses) can be added via params.
# Shadow=full ensures that only required outputs will be saved. 
rule get_virus_taxids:
    output: 
        "output/blast/10239.taxids"
    params:
       taxid = 10239
    conda:
        "https://raw.githubusercontent.com/avilab/virome-wrappers/master/blast/query/environment.yaml"
    resources:
        runtime = 120
    shell:
       "get_species_taxids.sh -t {params.taxid} > {output}"


rule megablast_virus:
    input:
        query = "output/{run}/repmaskedgood_{n}.fa",
        taxidlist = "output/blast/10239.taxids"
    output:
        out = temp("output/{run}/megablast-virus_{n}.tsv")
    params:
        program = "blastn",
        task = "megablast",
        db = "nt_v5",
        word_size = 16,
        evalue = 1e-6,
        outfmt = "'6 qseqid sacc staxid pident length qstart qend sstart send evalue'"
    threads: 8
    resources:
        runtime = lambda wildcards, attempt: attempt * 120,
        mem_mb = 26000
    wrapper:
        BLAST_QUERY

# Filter blastn hits for the cutoff value.
rule parse_megablast_virus:
    input:
        query = rules.megablast_virus.input.query,
        blast_result = rules.megablast_virus.output.out
    output:
        mapped = temp("output/{run}/megablast-virus_{n}_mapped.tsv"),
        unmapped = temp("output/{run}/megablast-virus_{n}_unmapped.fa")
    params:
        e_cutoff = 1e-6,
        outfmt = rules.megablast_virus.params.outfmt
    resources:
        runtime = lambda wildcards, attempt: attempt * 120,
        mem_mb = 4000
    wrapper:
        PARSE_BLAST


# Blastn, megablast and blastx input, output, and params keys must match commandline blast option names. 
# Please see https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a for all available options.
# Blast against nt virus database.
rule blastn_virus:
    input:
        query = rules.parse_megablast_virus.output.unmapped,
        taxidlist = "output/blast/10239.taxids"
    output:
        out = temp("output/{run}/blastn-virus_{n}.tsv")
    params:
        program = "blastn",
        db = "nt_v5",
        word_size = 11,
        evalue = 1e-6,
        outfmt = rules.megablast_virus.params.outfmt
    threads: 4
    resources:
        runtime = lambda wildcards, attempt: 840 + (attempt * 120),
        mem_mb = 26000
    wrapper:
        BLAST_QUERY

# Filter blastn hits for the cutoff value.
rule parse_blastn_virus:
    input:
        query = rules.blastn_virus.input.query,
        blast_result = rules.blastn_virus.output.out
    output:
        mapped = temp("output/{run}/blastn-virus_{n}_mapped.tsv"),
        unmapped = temp("output/{run}/blastn-virus_{n}_unmapped.fa")
    params:
        e_cutoff = 1e-6,
        outfmt = rules.megablast_virus.params.outfmt
    resources:
        runtime = lambda wildcards, attempt: attempt * 120,
        mem_mb = 4000
    wrapper:
        PARSE_BLAST

rule megablast_nt:
    input:
        query = rules.parse_blastn_virus.output.unmapped
    output:
        out = temp("output/{run}/megablast-nt_{n}.tsv")
    params:
        program = "blastn",
        task = "megablast",
        db = "nt_v5",
        word_size = 16,
        evalue = 1e-6,
        outfmt = rules.megablast_virus.params.outfmt
    threads: 8
    resources:
        runtime = lambda wildcards, attempt: attempt * 120,
        mem_mb = 30000
    wrapper:
        BLAST_QUERY

# Filter blastn hits for the cutoff value.
rule parse_megablast_nt:
    input:
        query = rules.megablast_nt.input.query,
        blast_result = rules.megablast_nt.output.out
    output:
        mapped = temp("output/{run}/megablast-nt_{n}_mapped.tsv"),
        unmapped = temp("output/{run}/megablast-nt_{n}_unmapped.fa")
    params:
        e_cutoff = 1e-6,
        outfmt = rules.megablast_virus.params.outfmt
    resources:
        runtime = lambda wildcards, attempt: attempt * 120,
        mem_mb = 4000
    wrapper:
        PARSE_BLAST


rule blastn_nt:
    input:
        query = rules.parse_megablast_nt.output.unmapped
    output:
        out = temp("output/{run}/blastn-nt_{n}.tsv")
    params:
        program = "blastn",
        db = "nt_v5",
        word_size = 11,
        evalue = 1e-6,
        outfmt = rules.megablast_virus.params.outfmt
    threads: 4
    resources:
        runtime = 400,
        mem_mb = 60000
    wrapper:
        BLAST_QUERY

# Filter blastn hits for the cutoff value.
rule parse_blastn_nt:
    input:
        query = rules.blastn_nt.input.query,
        blast_result = rules.blastn_nt.output.out
    output:
        mapped = temp("output/{run}/blastn-nt_{n}_mapped.tsv"),
        unmapped = temp("output/{run}/blastn-nt_{n}_unmapped.fa")
    params:
        e_cutoff = 1e-6,
        outfmt = rules.megablast_virus.params.outfmt
    resources:
        runtime = lambda wildcards, attempt: attempt * 120,
        mem_mb = 4000
    wrapper:
        PARSE_BLAST


# Blastx unmapped reads against nr virus database.
rule blastx_virus:
    input:
        query = rules.parse_blastn_nt.output.unmapped,
        taxidlist = "output/blast/10239.taxids"
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
        expand("output/{{run}}/{blastresult}_{{n}}_mapped.tsv", blastresult = BLAST)
    output:
        temp("output/{run}/all_{n}.csv")
    params:
        pp_sway = 0.5, 
        ranks_of_interest = RANKS_OF_INTEREST,
        dbfile = TAXON_DB
    resources:
        runtime = lambda wildcards, attempt: attempt * 120,
        mem_mb = 8000
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
    resources:
        runtime = lambda wildcards, attempt: attempt * 120,
        mem_mb = 8000
    run:
        tab = concatenate_tables(input, sep = ",", cols_to_integer = params.ranks)
        mask = tab.superkingdom == VIRUSES_TAXID
        mask = mask.fillna(False)
        vir = tab[mask]
        non_vir = tab[~mask]
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
