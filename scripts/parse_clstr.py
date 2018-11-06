
from re import match, findall, sub, compile
import scipy.stats as ss

path = "/Users/taavi/Downloads/SRR5580355_cdhit.fa.clstr_short"

@do_cprofile
def iterate_clstr(clstr_path, clstr_ids, top_n = 3):
  reid = compile(r"(?<=[>])\w+.\d+")
  resim = compile(r"\*$|\d+\.\d+(?=%)")
  with open(clstr_path, "r") as clstr_file, open(clstr_ids, "w") as ids_out:
    for line in clstr_file:
      if match(r"(>Clus)", line):
        try:
          cls
        except NameError:
          cls = []
        else:
          ids = [reid.findall(record)[0] for record in cls]
          similarity = [resim.findall(record)[0] for record in cls]
          ids_filtered = [tup for tup in zip(ids, similarity) if tup[1] not in set(["*", "100.00"])]
          ranks = list(ss.rankdata([tup[1] for tup in ids_filtered], method = "dense"))
          parsed_ids = [tup[0][0] for tup in zip(ids_filtered, ranks) if tup[1] <= top_n]
          ids_out.writelines("%s\n" % k for k in parsed_ids)
          cls = []
      if match(r"\d+", line):
        cls.append(line)


iterate_clstr(path, "/Users/taavi/Downloads/SRR5580355_cdhit.ids")

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
result = iterate_clstr(path, "/Users/taavi/Downloads/SRR5580355_cdhit.ids")
