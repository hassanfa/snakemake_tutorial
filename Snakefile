# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

result_dir = "result"
fastq_dir = "fastq"
analysis_dir = os.path.join(result_dir, "analysis")
ref_dir = os.path.join(result_dir, "ref")
log_dir = os.path.join(result_dir, "log")
benchmark_dir = os.path.join(result_dir, "benchmark")

genome_ref = "ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.22.fa.gz"

rule all:
  input:
    
#    multiext(os.path.join(ref_dir, "Homo_sapiens.GRCh37.75.dna.chromosome.22.fa"), ".fai", ".amb", ".ann", ".bwt", ".pac", ".sa") 

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

rule samtools_index:
  input:
    os.path.join(ref_dir, "Homo_sapiens.GRCh37.75.dna.chromosome.22.fa")
  output:
    os.path.join(ref_dir, "Homo_sapiens.GRCh37.75.dna.chromosome.22.fa.fai")
  version: "0.0.1"
  message: "Preparing samtools index for reference"
  log: os.path.join(log_dir, "reference_index.log")
  benchmark: os.path.join(benchmark_dir, "reference_index.tsv")
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
  shell:
    """
bwa index {input} &> {log}
    """ 

rule mapping:
  input:
    fa = rules.download_ref.output,
    fa_index = rules.bwa_index.output,
    reads = expand(os.path.join(fastq_dir, "{sample}_{rp}.fastq", rp=["1", "2"]))
  output:
    os.path.join(analysis_dir, "{sample}.sorted.bam")
  version: "0.0.1"
  message: "aligning to genome"
  log: os.path.join(log_dir, "mapping_{sample}.log")
  benchmark: os.path.join(benchmark_dir, "mapping_{sample}.tsv")
  shell:
    """
bwa mem {input.fa} {input.reads} | samtools sort | samtools view -Sb - > {output}
    """
