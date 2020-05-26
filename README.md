#### Requirements

Make sure you have conda installed.

#### Step 1: Create a conda environment

```
conda create -n cg_development_day -c bioconda -c conda-forge python=3.6 pip
```

#### Step 2: Activate conda environment

```
conda activate cg_development_day
```

#### Step 3: Install all Snakemake dependencies via pip

```
pip install -r requirements.txt 
```

### Conda stuff 
#### Step 1: Install required softwares 

```
conda install -c bioconda -c conda-forge vardict-java bcftools samtools bwa rtg-tools
```


#### Step 2: Run workflow

```
snakemake --cores 2 -p
```

#### Step 3: Generate a report

```
snakemake --report && open report.html
```
