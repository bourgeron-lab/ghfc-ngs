# Using **ghfc-ngs** from Institut Pasteur

This workflow was developed to be used in the Institut Pasteur cluster. Here are the main considerations specific to this cluster.

## Installation

To avoid having to maintain the cluster version up-to-date, nextflow allow to run a workflow directly from a github repository. He we will also extend the principle to the running script of the pipeline, thus only needing a small script that can be copy-pasted from here and that should not need to be updated later on.

Create a script called `ghfc-ngs`, make it executable (`chmod +x`) and store it somewhere in your path. Copy the following content in this script.

```bash
#!/bin/bash

bash -c "$(curl -s https://raw.githubusercontent.com/bourgeron-lab/ghfc-ngs/refs/heads/main/run_pipeline.sh )" -s $@
```

By default, this script exists in the `/pasteur/helix/projects/ghfc_wgs/tools/bin` directory.

If not already done, you can add this directory to your PATH by adding the following line in your `~/.bashrc` file:

```bash
PATH="/pasteur/helix/projects/ghfc_wgs/tools/bin:$PATH"
```

## Usage

The ghfc-ngs workflow requires at least a parameter file, and has several optional parameters.

If you want to display the usage message:
`ghfc-ngs -h`

To create a parameter files, the easiest way is to copy the one from this repository ([available here](params.yml)) and then edit it with your favorite code editor.

```bash
curl xxxxx
```

## Tuning

Several parameters could be tuned, but the default values should be good for Pasteur's cluster.

### 1. Nextflow config file (*nextflow.config*)

> **Important**
>
> the nextflow.config from the github is used by default, but you can override it with the `--config` option of the runner.

In this file, the user can define the default resources for a process (*e.g.* currently 1 cpu and 4 GB of memory). It is also possible to adjust process-specific requirements.

```bash
withName: 'BWA_MEM2_ALIGN' {
        cpus = 95
        memory = '460.GB'
    }
```

Additionaly the container images used by the workflow are defined there and can be changed or updated (for newer version, compatibility issues or unavaibility of one).

```bash
# DeepVariant process settings requires version > 1.8
withName: 'DV_MAKE_EXAMPLES' {
        cpus = 95
        memory = '460.GB'
        container = 'docker://google/deepvariant:1.9.0'
    }
```

If there is a minimal version or a broking change identified in a different version, it is a good idea to indicate it as an issue in the github.

The slurm parameters (such as the partition or the qos) are also defined in this config file.

```bash
slurm {
        process.executor = 'slurm'
        process.clusterOptions = '-p ghfc --qos=ghfc'
    }
```

### 2. Workflow parameter file (*params.yml*)

The tuning in this file is mostly restricted to the maximum resources available per node.

```bash
# Resource limits
max_memory: "460.GB"
max_cpus: 95
max_time: "240.h"
```

### 3. hardcoded in the run script (*run_pipeline.sh*)

A few interesting things are hardcoded in the run script, but they are specific to Pasteur cluster and most likely do not need to be adjusted.

There is the loading of the required modules of the cluster.

```bash
module load graalvm
module load apptainer
module load graphviz
module load nextflow
```

The ulimit parameters and OPENBLAS_NUM_THREADS are here to solve numerous problems mainly in DeepVariant with tensorflow, or for GLNexus for its RockDB use.

```bash
ulimit -v unlimited
ulimit -Sn 65536
ulimit -u 65536
export OPENBLAS_NUM_THREADS=1
```
