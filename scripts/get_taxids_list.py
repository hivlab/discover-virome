from snakemake.shell import shell

taxdict = snakemake.params

for k,v in taxdict.items():
    shell("get_species_taxids.sh -t {} > blast/{}.taxid".format(v, k))
    neg_taxid_files = ["blast/{}.taxid".format(k) for k,v in taxdict.items() if k != "viruses"]

shell("cat {neg_taxid_files} > blast/negative.taxid")
