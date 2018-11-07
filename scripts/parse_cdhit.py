
from re import match, findall, compile
import scipy.stats as ss
from Bio import SeqIO
from common.helpers import subset_records

clstr = "/Users/taavi/Downloads/SRR5580355_cdhit.fa.clstr_short"
topn_clstr = "/Users/taavi/Downloads/SRR5580355_cdhit.ids"
args = list([clstr, topn_clstr])

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
records = subset_records(SeqIO.parse(snakemake.input[2], "fasta"), set(topn_ids))

# Write topn subset to file
count = SeqIO.write(records, snakemake.output[1], 'fasta')

