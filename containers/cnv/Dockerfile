# What the cnv_coverage and cnv_panel steps run with: mosdepth for the binned depth of each CRAM,
# and python3 with numpy, pandas and pyBigWig for modules/cnv/scripts/cnv_coverage, which turns it
# into the normalised bigWigs, builds the control panel and fits a cohort against it.
#
# mosdepth is pinned to the version of the MOSDEPTH process's biocontainer, so that the two
# coverage files of a sample are made by the same tool.
#
#   docker build --platform linux/amd64 -t fcliquet/cnv-coverage:1.0.0 containers/cnv
#   docker push fcliquet/cnv-coverage:1.0.0
#
# procps is there for Nextflow, which reads task metrics with ps.
FROM mambaorg/micromamba:2.0.5

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        "mosdepth=0.3.11" \
        "python>=3.12,<3.13" \
        "numpy>=2,<3" \
        "pandas>=2.2,<3" \
        "pybigwig>=0.3.22,<0.4" \
        procps-ng \
    && micromamba clean -a -y

# Apptainer does not run the micromamba entrypoint, so the environment must be on PATH as is
ENV PATH=/opt/conda/bin:${PATH}

# Fail the build rather than the first task if anything did not resolve, or if pyBigWig was
# built without numpy, which cnv_coverage writes and reads its bins through
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
RUN mosdepth --version && python3 --version && ps --version \
    && python3 -c 'import numpy, pandas, pyBigWig; assert pyBigWig.numpy; print(numpy.__version__, pandas.__version__)'
