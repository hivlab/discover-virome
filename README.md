
[![Travis-CI Build Status](https://travis-ci.org/<USERNAME>/<REPO>.svg?branch=master)](https://travis-ci.org/<USERNAME>/<REPO>)

The goal of this repo is to more reproducibly recreate [VirusSeeker Virome workflow](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5326578/). Original code is available on github.com/guoyanzhao.


## Setup environment and install prerequisites

### Download RepBase for RepeatMasker

Obtain access to RepBase from www.girinst.org. 
Download RepBase into separate directory e.g. "databases/repbase" under your home directory or whatever other location is good for you. 
For downloading set environment variables for GIRUSER and GIRPASS. 
```
mkdir databases/repbase
cd databases/repbase

GIRUSER=<your-gir-user-name>
GIRPASS=<your-gir-password>

wget --user $GIRUSER --password $GIRPASS --no-check-certificate http://www.girinst.org/server/RepBase/protected/repeatmaskerlibraries/RepBaseRepeatMaskerEdition-20170127.tar.gz

gunzip RepBaseRepeatMaskerEdition-20170127.tar.gz
tar xvf RepBaseRepeatMaskerEdition-20170127.tar
rm RepBaseRepeatMaskerEdition-20170127.tar
```

### Install miniconda

Download and install miniconda https://conda.io/docs/user-guide/install/index.html.
In case of Linux, following should work:
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

### Install environment

Create conda environment with preinstalled **snakemake**:
```
conda create -n snakemake -c bioconda -c conda-forge snakemake
```

### Activate environment

```
source activate snakemake
```

### Clone this repo and cd to repo
(Change URL accordingly if using HTTPS)

```
git clone git@github.com:avilab/vs.git
```

## 


## Example

### Dry run

```
snakemake -n
```

### Create workflow graph

```
snakemake --dag | dot -Tsvg > graph/dag.svg
```

### Run workflow

This workflow is designed to be run in cluster. `cluster.json` configuration file may need some customisation, for example partition name. Memory nad maximum runtime limits are optimised for 100 splits. Number of splits can be specified in `config.yaml` file with n_files option (currently n_files is 2). Installation of software dependencies is taken care by conda, hence there is software installation overhead when you run this workflow for the first time in new environment. 

Example workflow submission script for slurm cluster, where values for job name, cluster partition name, time and memory constraints, and slurm log path (output) are taken from cluster.json: 
```
snakemake -j --use-conda --cluster-config cluster.json  \
             --cluster "sbatch -J {cluster.name} \
             -p {cluster.partition} \
             -t {cluster.time} \
             --mem {cluster.mem} \
             --output {cluster.output}"
```

You may want to use also following flags when running this workflow in cluster:
```
--max-jobs-per-second 1 --max-status-checks-per-second 10 --rerun-incomplete --keep-going
```

All other possible [snakemake execution](https://snakemake.readthedocs.io/en/stable/executable.html) options can be printed by calling `snakemake -h`.

### Exit/deactivate environment

Conda environment can be closed with the following command when work is finished:
```
source deactivate
```

