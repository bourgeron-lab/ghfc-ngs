# What the melt step's processes run with, MELT itself excepted: Java 11 to run MELT.jar,
# samtools and bowtie2 that MELT calls out to, bcftools/htslib to merge the family VCFs and to
# read gnomAD-SV, and python3 for the cohort scoring.
#
# MELT is not in this image, and must never be: its licence (University of Maryland, academic
# use) is non-exclusive and non-transferable, so it cannot be redistributed in an image anyone
# can pull. The MELT distribution is staged into each task from params.melt_dir instead.
#
# The versions are those the cohort's first 4,018 calls were made with on Maestro: OpenJDK 11,
# samtools 1.21, bowtie2 2.5.4.
#
#   docker build --platform linux/amd64 -t fcliquet/melt-deps:1.0.0 containers/melt
#   docker push fcliquet/melt-deps:1.0.0
#
# procps is there for Nextflow, which reads task metrics with ps.
FROM mambaorg/micromamba:2.0.5

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        "openjdk>=11,<12" \
        "samtools=1.21" \
        "bowtie2=2.5.4" \
        "htslib=1.21" \
        "bcftools=1.21" \
        "python>=3.11,<3.13" \
        procps-ng \
    && micromamba clean -a -y

# Apptainer does not run the micromamba entrypoint, so the environment must be on PATH as is
ENV PATH=/opt/conda/bin:${PATH}

# Fail the build rather than the first task if anything did not resolve. pipefail, or a missing
# samtools would pass through `| head -1`
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
RUN java -version 2>&1 | head -1 && samtools --version | head -1 && bowtie2 --version | head -1 \
    && bcftools --version | head -1 && python3 --version && ps --version
