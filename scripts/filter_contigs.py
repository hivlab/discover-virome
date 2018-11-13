
from Bio import SeqIO
import re

cov_pattern = re.compile("cov_([0-9.]+)")
min_length, min_coverage = snakemake.params

with open(snakemake.input[0], 'rU') as input_fasta, open(snakemake.output[0], 'w') as filtered_fasta:
  for contig in SeqIO.parse(input_fasta, 'fasta'):
    coverage = cov_pattern.search(contig.name)
		if len(contig) >= min_length and float(result.group(1)) >= min_coverage:
		  SeqIO.write(contig, filtered_fasta, 'fasta')
