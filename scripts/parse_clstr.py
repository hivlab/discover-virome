
from re import match, findall, sub, compile
import scipy.stats as ss

path = "/Users/taavi/Downloads/SRR5580355_cdhit.fa.clstr"

@do_cprofile
def iterate_clstr(clstr_path, clstr_ids, top_n = 3):
  combo = compile(r"(?<=[>])\w+.\d+|\*$|\d+\.\d+(?=%)")
  clus = compile(r"(>Clus)")
  with open(clstr_path, "r") as clstr_file, open(clstr_ids, "w") as ids_out:
    for line in clstr_file:
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


iterate_clstr(path, "/Users/taavi/Downloads/SRR5580355_cdhit.ids")

reid = compile(r"(?<=[>])\w+.\d+")
resim = compile(r"\*$|\d+\.\d+(?=%)")
combo = compile(r"(?<=[>])\w+.\d+|\*$|\d+\.\d+(?=%)")
cls = ['141nt, >SRR5580355.436734... at 141:1:34:174/-/99.29%', '141nt, >SRR5580355.171785... at 141:1:78:218/-/99.29%']

m = [tuple(combo.findall(record)) for record in cls]


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
