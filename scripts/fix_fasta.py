from Bio import SeqIO

fixed_seqs = []
for record in SeqIO.parse(snakemake.input[0], "fasta"):
    record.description = ""
    fixed_seqs.append(record)

SeqIO.write(fixed_seqs, snakemake.output[0], "fasta")
