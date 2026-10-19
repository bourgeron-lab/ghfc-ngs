# What the sampleqc step runs with: somalier to extract each CRAM's allele counts at its sites
# and to relate the samples, bcftools to stream a family's norm.bcf into the Mendelian counts, and
# python3, standard library only, for modules/sampleqc/scripts/sampleqc_classify.
#
#   docker build --platform linux/amd64 -t fcliquet/sampleqc:1.0.0 containers/sampleqc
#   docker push fcliquet/sampleqc:1.0.0
#
# procps is there for Nextflow, which reads task metrics with ps.
FROM mambaorg/micromamba:2.0.5

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        "somalier=0.3.5" \
        "htslib=1.24" \
        "bcftools=1.24" \
        "python>=3.12,<3.13" \
        procps-ng \
    && micromamba clean -a -y

# Apptainer does not run the micromamba entrypoint, so the environment must be on PATH as is
ENV PATH=/opt/conda/bin:${PATH}

# Fail the build rather than the first task if anything did not resolve. somalier prints its
# version on stderr and exits 0 with no subcommand.
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
RUN somalier 2>&1 | grep -q 'somalier version: 0.3.5' && bcftools --version | head -1 \
    && python3 --version && ps --version && base64 --version | head -1
