from snakemake.shell import shell
from os.path import dirname 

outdir = dirname(snakemake.output[0])
if len(outdir) == 0:
    outdir = "."
taxdict = snakemake.params

for k,v in taxdict.items():
    shell("get_species_taxids.sh -t {} > {}/{}.taxid".format(v, outdir, k))

neg_taxid_files = ["{}/{}.taxid".format(outdir, k) for k,v in taxdict.items() if k != "viruses"]
shell("cat {neg_taxid_files} > {outdir}/negative.taxid")
