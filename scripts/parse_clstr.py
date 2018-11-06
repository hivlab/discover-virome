
from re import search, findall, sub, compile
import scipy.stats as ss

path = "/Users/taavi/Downloads/SRR5580355_cdhit.fa.clstr_short"

@do_cprofile
def iterate_clstr(path, top_n = 3):
  reid = compile(r"(?<=[>])\w+.\d+")
  resim = compile(r"\*$|\d+\.\d+(?=%)")
  parsed = []
  with open(path) as handle:
    for line in handle:
      if search(r"(^>Clus)", line):
        try:
          cls
        except NameError:
          cls = []
        else:
          ids = [reid.findall(record)[0] for record in cls]
          similarity = [resim.findall(record)[0] for record in cls]
          ids_filtered = [tup for tup in zip(ids, similarity) if tup[1] not in set(["*", "100.00"])]
          ranks = list(ss.rankdata([tup[1] for tup in ids_filtered], method = "dense"))
          parsed.append([tup[0][0] for tup in zip(ids_filtered, ranks) if tup[1] <= top_n])
          cls = []
      if search(r"^\d+", line):
        cls.append(line)
  return sum(parsed, [])


out = iterate_clstr(path)

import cProfile

def do_cprofile(func):
    def profiled_func(*args, **kwargs):
        profile = cProfile.Profile()
        try:
            profile.enable()
            result = func(*args, **kwargs)
            profile.disable()
            return result
        finally:
            profile.print_stats()
    return profiled_func

# perform profiling
result = iterate_clstr(path)
