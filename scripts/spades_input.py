
from Bio import SeqIO
import gzip
from common.helpers import get_ids

# Save joined- and unjoined sequences into separate files for spades
joined = SeqIO.parse(gzip.open(snakemake.input[0], "rt"), "fastq")
joined_ids = get_ids(joined)

with open(snakemake.output[0], "w") as join, open(snakemake.output[1], "w") as un:
  for record in SeqIO.parse(snakemake.input[1], "fastq"):
        if record.id in joined_ids:
            SeqIO.write(record, join, 'fastq')
        else:
            SeqIO.write(record, un, 'fastq')
