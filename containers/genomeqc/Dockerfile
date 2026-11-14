# What the genomeqc step runs with: samtools to count a CRAM's reads on a few contigs and take its
# stats on chr20, bcftools to stream a family's norm.bcf, and python3, standard library only, for
# modules/genomeqc/scripts/genomeqc_tables. htslib 1.24, as in fcliquet/sampleqc.
#
#   docker build --platform linux/amd64 -t fcliquet/genomeqc:1.0.0 containers/genomeqc
#   docker push fcliquet/genomeqc:1.0.0
#
# procps is there for Nextflow, which reads task metrics with ps.
FROM mambaorg/micromamba:2.0.5

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        "samtools=1.24" \
        "htslib=1.24" \
        "bcftools=1.24" \
        "python>=3.12,<3.13" \
        procps-ng \
    && micromamba clean -a -y

# Apptainer does not run the micromamba entrypoint, so the environment must be on PATH as is
ENV PATH=/opt/conda/bin:${PATH}

# Fail the build rather than the first task if anything did not resolve
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
RUN samtools --version | head -1 | grep -q 'samtools 1.24' && bcftools --version | head -1 \
    && python3 --version && ps --version && base64 --version | head -1 && awk 'BEGIN { exit 0 }'
