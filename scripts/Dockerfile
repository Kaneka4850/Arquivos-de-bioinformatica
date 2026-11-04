FROM mambaorg/micromamba:1.5.10
RUN micromamba install -y -n base -c conda-forge -c bioconda \
      python=3.11 samtools bcftools htslib bwa art && micromamba clean -a -y
WORKDIR /data
