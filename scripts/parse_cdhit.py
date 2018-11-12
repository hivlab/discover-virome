
from re import match, findall, compile
import scipy.stats as ss
from Bio import SeqIO
from subprocess import Popen, PIPE
import gzip
from common.helpers import subset_records
from common.helpers import get_ids

def iterate_clstr(clstr, topn_clstr, top_n = 3):
  combo = compile(r"(?<=[>])\w+.\d+|\*$|\d+\.\d+(?=%)")
  clus = compile(r"(>Clus)")
  with open(clstr, "r") as clusters, open(topn_clstr, "w") as ids_out:
    for line in clusters:
      if clus.match(line):
        try:
          cls
        except NameError:
          cls = []
        else:
          matches = [tuple(combo.findall(record)) for record in cls]
          ids_filtered = [tup for tup in matches if tup[1] not in ["*", "100.00"]]
          ranks = list(ss.rankdata([tup[1] for tup in ids_filtered], method = "dense"))
          parsed_ids = [tup[0][0] for tup in zip(ids_filtered, ranks) if tup[1] <= top_n]
          ids_out.writelines("%s\n" % k for k in parsed_ids)
          cls = []
      else:
        cls += [line]

# Parse and save cluster ids
iterate_clstr(snakemake.input[1], snakemake.output[0], snakemake.params[0])

# Import parsed cluster ids
with open(snakemake.output[0]) as seq_ids:
  topn_ids = seq_ids.read().splitlines()

# Subset sequences using topn ids
records = subset_records(SeqIO.parse(snakemake.input[2], "fastq"), set(topn_ids))

# Write topn subset to file
count = SeqIO.write(records, snakemake.output[1], 'fasta')

# Concatenate representative sequences and topn
cmd = "cat {} {}".format(snakemake.input[0], snakemake.output[1])
with open(snakemake.output[2], "w") as merged:
  process = Popen(cmd.split(' '), stdout = merged, stderr = PIPE)
  stdout, stderr = process.communicate()

# Save joined- and unjoined sequences into separate files for spades
joined = SeqIO.parse(gzip.open(snakemake.input[3], "rt"), "fastq")
joined_ids = get_ids(joined)

with open(snakemake.output[3]) as joned, open(snakemake.output[4]) as un:
  for record in SeqIO.parse(snakemake.output[2], "fasta"):
        if record.id in joined_ids:
            SeqIO.write(record, joined, 'fasta')
        else:
            SeqIO.write(record, un, 'fasta')
