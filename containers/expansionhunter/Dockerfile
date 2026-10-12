# ExpansionHunter v5.0.0 for the expansionhunter step, with what its three processes need
# around it: samtools to sort and index the realigned BAM, bcftools/htslib to merge the family
# and cohort VCFs, and python3 for the genotype tables.
#
# Bioconda build 5 is pinned on purpose. It links htslib 1.22, which reads CRAM 3.1; Illumina's
# release binary embeds an htslib older than 1.14 and fails on CRAM 3.1 with "Failure to decode
# slice", which is why the cohort's first run had to transcode every CRAM 3.1 to 3.0 first.
#
#   docker build --platform linux/amd64 -t fcliquet/expansionhunter:5.0.0 containers/expansionhunter
#   docker push fcliquet/expansionhunter:5.0.0
#
# procps is there for Nextflow, which reads task metrics with ps.
FROM mambaorg/micromamba:2.0.5

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        expansionhunter=5.0.0=hc26b3af_5 \
        "htslib>=1.22.1,<1.23" \
        "samtools>=1.22,<1.23" \
        "bcftools>=1.22,<1.23" \
        "python>=3.11,<3.13" \
        procps-ng \
    && micromamba clean -a -y

# Apptainer does not run the micromamba entrypoint, so the environment must be on PATH as is
ENV PATH=/opt/conda/bin:${PATH}

# Fail the build rather than the first task if anything did not resolve. pipefail, or a missing
# samtools would pass through `| head -1`
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
RUN ExpansionHunter --version && samtools --version | head -1 && bcftools --version | head -1 \
    && python3 --version && ps --version
