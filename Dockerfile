FROM hassanf/miniconda3_4.6.14:latest

LABEL maintainer="Hassan Foroughi hassan dot foroughi at scilifelab dot se" 
LABEL description="Bioinfo tools for snakemake tutorial"
LABEL version="0.0.1"

# create necessary directories
# install balsamic and it's environments
# symlink libreadline for picard to function properly 
RUN yum install -y which wget bzip2 graphviz graphviz-devel git gcc fontconfig a& \
    export PATH=/usr/local/miniconda/bin:$PATH && \
    export LC_ALL=en_US.utf-8 && \
    export LANG=en_US.utf-8 && \
    conda install -c bioconda -c conda-forge vardict-java bcftools samtools bwa rtg-tools python=3.6 pip && \ 
    pip install -r https://raw.githubusercontent.com/hassanfa/snakemake_tutorial/master/requirements.txt && \
    conda clean --all -y

# The following fixes the error for Click
# RuntimeError: Click will abort further execution because Python 3 was
# configured to use ASCII as encoding for the environment. Consult
# https://click.palletsprojects.com/en/7.x/python3/ for mitigation steps.
ENV LC_ALL=en_US.utf-8
ENV LANG=en_US.utf-8
ENV PATH="/usr/local/miniconda/bin:${PATH}"
