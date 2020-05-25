# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

result_dir = "result"
fastq_dir = "fastq"
analysis_dir = os.path.join(result_dir, "analysis")
ref_dir = os.path.join(result_dir, "ref")
sdf_dir = os.path.join(ref_dir, "sdf")
log_dir = os.path.join(result_dir, "log")
benchmark_dir = os.path.join(result_dir, "benchmark")

genome_ref = "ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.22.fa.gz"
ref_region = "resources/region.bed"
vcf_truth = "resources/truth_set.vcf.gz"
samples = ["tumor", "normal"]
VAF_threshold = 0.01

rule all:
  input:
    os.path.join(analysis_dir, "somatic.vcf.gz"),
    os.path.join(analysis_dir, "var_eval_af", "summary_af.svg"),
    os.path.join(analysis_dir, "var_eval_dp", "summary_dp.svg"),
    expand(os.path.join(analysis_dir, "{sample}.sorted.bam"), sample = samples)

rule download_ref:
  output:
    os.path.join(ref_dir, "Homo_sapiens.GRCh37.75.dna.chromosome.22.fa")
  params: genome_ref_path = genome_ref
  version: "0.0.1" 
  message: "Downloading reference data"
  log: os.path.join(log_dir, "download_ref.log")
  benchmark: os.path.join(benchmark_dir, "download_ref.tsv")
  shell:
    """
wget -a {log} -O - {params.genome_ref_path} | gunzip > {output} 
    """

rule ref_format:
  input:
    os.path.join(ref_dir, "Homo_sapiens.GRCh37.75.dna.chromosome.22.fa")
  output:
    directory(os.path.join(ref_dir, sdf_dir))
  version: "0.0.1"
  message: "Generating SDF file for reference Fasta"
  log: os.path.join(log_dir, "ref_format.log")
  benchmark: os.path.join(benchmark_dir, "ref_format.tsv")
  conda: "envs/rtg.yaml"
  shell:
    """
rtg format -o {output} {input}
    """

rule samtools_index:
  input:
    os.path.join(ref_dir, "Homo_sapiens.GRCh37.75.dna.chromosome.22.fa")
  output:
    os.path.join(ref_dir, "Homo_sapiens.GRCh37.75.dna.chromosome.22.fa.fai")
  version: "0.0.1"
  message: "Preparing samtools index for reference"
  log: os.path.join(log_dir, "reference_index.log")
  benchmark: os.path.join(benchmark_dir, "reference_index.tsv")
  conda: "envs/bwa.yaml"
  shell:
    """
samtools faidx {input} &> {log} 
    """ 

rule bwa_index:
  input:
    os.path.join(ref_dir, "Homo_sapiens.GRCh37.75.dna.chromosome.22.fa")
  output:
    multiext(os.path.join(ref_dir, "Homo_sapiens.GRCh37.75.dna.chromosome.22.fa"), ".amb", ".ann", ".bwt", ".pac", ".sa")
  version: "0.0.1"
  message: "Preparing bwa index for reference"
  log: os.path.join(log_dir, "reference_index.log")
  benchmark: os.path.join(benchmark_dir, "reference_index.tsv")
  singularity: "docker://clinicalgenomics/bwa:0.7.17"
  conda: "envs/bwa.yaml"
  shell:
    """
bwa index {input} &> {log}
    """ 

def get_fastq_files(wildcards):
  """Return a list of fastq files for a sample"""
  return expand(os.path.join(fastq_dir, "{sample}_{readpair}.fastq"), readpair=[1, 2], **wildcards) 

rule mapping:
  input:
    fa = rules.download_ref.output,
    fa_index = rules.bwa_index.output,
    reads = get_fastq_files
  output:
    os.path.join(analysis_dir, "{sample}.sorted.bam")
  version: "0.0.1"
  message: "aligning to genome"
  log: os.path.join(log_dir, "mapping_{sample}.log")
  benchmark: os.path.join(benchmark_dir, "mapping_{sample}.tsv")
  conda: "envs/bwa.yaml"
  shell:
    """
(bwa mem {input.fa} {input.reads} | samtools sort | samtools view -Sb - > {output} && samtools index {output} ) &> {log}
    """

rule vardict_call:
  input:
    fa=rules.download_ref.output,
    fai=rules.samtools_index.output,
    tumorbam=os.path.join(analysis_dir, "tumor.sorted.bam"),
    normalbam=os.path.join(analysis_dir, "normal.sorted.bam"),
    bed= ref_region
  output:
    vcf=os.path.join(analysis_dir, "somatic.vcf.gz")
  params:
    tumorname="tumor",
    normalname="normal",
    VAF=VAF_threshold,
    vardict='-c 1 -S 2 -E 3 -g 4'
  version: "0.0.1"
  message: "Variant calling"
  log: os.path.join(log_dir, "vardict_call.log")
  benchmark: os.path.join(benchmark_dir, "vardict_call.tsv")
  conda: "envs/vardict.yaml"
  shell:
    """
vardict-java \
  -G {input.fa} \
  -f {params.VAF} \
  -N tumor \
  -b \"{input.tumorbam}|{input.normalbam}\" \
  {params.vardict} \
  {input.bed} \
  | testsomatic.R \
  | var2vcf_paired.pl -N \"{params.tumorname}|{params.normalname}\" \
    -f {params.VAF} \
  | bcftools filter -e \'STATUS !~ \".*Somatic\"\' \
  | bcftools view -f PASS --output-file {output.vcf} --output-type z;
tabix -p vcf {output.vcf};
    """ 

rule eval_variantcall:
  input:
    call = os.path.join(analysis_dir, "somatic.vcf.gz"),
    truth = vcf_truth, 
    sdf = rules.ref_format.output
  output:
    af_summary = report(os.path.join(analysis_dir, "var_eval_af", "summary_af.svg"), caption="resources/AF_eval.rst", category="Evaluation"),
    dp_summary = report(os.path.join(analysis_dir, "var_eval_dp", "summary_dp.svg"), caption="resources/DP_eval.rst", category="Evaluation"),
  params:
    af_eval = os.path.join(analysis_dir, "var_eval_af"),
    dp_eval = os.path.join(analysis_dir, "var_eval_dp"),
  version: "0.0.1"
  message: "Evaluating variant call"
  log: os.path.join(log_dir, "eval_variantcall.log")
  benchmark: os.path.join(benchmark_dir, "eval_variantcall.tsv")
  conda: "envs/rtg.yaml"
  shell:
    """
rm -r {params.af_eval} && rm -r {params.dp_eval};
rm {output.af_summary} && rm {output.dp_summary};
rtg vcfeval -b {input.truth} -c {input.call} \
  -o {params.af_eval} -t {input.sdf} \
  --sample SAMPLE,tumor --vcf-score-field INFO.AF &> {log};
rtg rocplot --svg {output.af_summary} {params.af_eval}/weighted_roc.tsv.gz &>> {log}; 
rtg vcfeval -b {input.truth} -c {input.call} \
  -o {params.dp_eval} -t {input.sdf} \
  --sample SAMPLE,tumor --vcf-score-field INFO.DP &>> {log}
rtg rocplot --svg {output.dp_summary} {params.dp_eval}/weighted_roc.tsv.gz &>> {log}; 
    """
