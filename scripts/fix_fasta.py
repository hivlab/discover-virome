from Bio import SeqIO

fixed_seqs = []
for record in SeqIO.parse(snakemake.input[0], "fasta"):
    record.description = ""
    record.id  = f"{snakemake.params[0]}_{record.id}"
    fixed_seqs.append(record)

SeqIO.write(fixed_seqs, snakemake.output[0], "fasta")
