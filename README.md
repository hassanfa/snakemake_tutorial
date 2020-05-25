#### Requirements

Make sure you have conda installed.

#### Step 1: Create a conda environment

```
conda create -n cg_development_day -c bioconda -c conda-forge python=3.6 pip
```

#### Step 3: Activate conda environment

```
conda activate cg_development_day
```

#### Step 4: Install all Snakemake dependencies via pip

```
pip install -r requirements.txt 
```

#### Step 5: Install required softwares 

```
conda install -c bioconda -c conda-forge vardict-java bcftools samtools bwa rtg-tools
```


#### Step 6: Run workflow

```
snakemake --cores 2 -p
```

#### Step 7: Generate a report

```
snakemake --report && open report.html
```
