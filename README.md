# negative-control-genomicc

### Motivation
This is a repository used to generate a biological negative control for the GenOMICC cohort. It takes as input BGEN files for a given cohort, an input list of binding QTLs (bQTLs) and trans-actors (ex: expression QTLs; eQTLs), as well as regions of the genome to exclude, and pick a randomly-sampled MAF-matched SNPs for the input SNPs you specify (bQTLs + eQTLs).

The idea behind this is to test the biological hypothesis that epistasis through altered binding of a transcription factor (TF) and a trans-acting QTL which regulates binding at a different level (for example, by expression of that TF) drives complex trait variation in a genomic cohort.

### Requirements
You must have nextflow and singularity installed. This can be achieved through conda and the `env.yaml` file provided in this repository.

#### Conda
You can install conda by the following the directions here: https://conda.io/projects/conda/en/latest/user-guide/install/index.html

Example:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

#### Nextflow & Singularity
Once installed, you can create the nextflow environment by:
```
conda env create --file env.yaml
```

### Setup
You must provide the following information in a run-specific configuration file (example: `genomicc.config`)

`OUTDIR`            : Path to output directory

`BQTLS`             : Path to CSV containing the bQTLs you would like to find negative control SNPs for

`EQTLS`             : Path to CSV containing the eQTLs you would like to find negative control SNPs for

`ASSEMBLY`          : Assembly name (either hg19, grch37, hg38 or grch38)

`EXCLUSION_REGIONS` : Path to CSV file containing the LD blocks for lead SNPs, which you would like not to sample from

`CONFIG_FILE`       : Config file for targene test run

`ESTIMANDS_FILE`    : Estimands file for targene test run

`BGEN_FILES`        : /path/to/BGEN/cohort_chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}.{bgen,bgen.bgi,sample}"

### Run
To run the pipeline, start a screen session, activate your `nextflow` environment and then run with the profile for your platform (in my case, Eddie):

```
screen
conda activate nextflow
nextflow run main.nf -profile eddie -resume
```

### Output
Results will be stores in the OUTDIR specified in your configuration above. 