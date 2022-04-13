# ![nf-core/circdna](docs/images/nf-core/circdna_logo_light.png#gh-light-mode-only) ![nf-core/circdna](docs/images/nf-core/circdna_logo_dark.png#gh-dark-mode-only)

[![GitHub Actions CI Status](https://github.com/nf-core/circdna/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/circdna/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/nf-core/circdna/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/circdna/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/circdna/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23circdna-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/circdna)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->

**nf-core/circdna** is a bioinformatics best-practice analysis pipeline for the identification of circular DNAs in eukaryotic cells. The pipeline is able to process WGS, ATAC-seq data or Circle-Seq data generated from short-read sequencing technologies.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/circdna/results).

## Pipeline summary

1. Merge re-sequenced FastQ files ([`cat`](http://www.linfo.org/cat.html))
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Adapter and quality trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
4. Map reads using BWA-MEM ([`BWA`](https://github.com/lh3/bwa))
5. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
6. Choice of multiple circular DNA identification routes
   1. [`Circle-Map ReadExtractor`](https://github.com/iprada/Circle-Map) -> [`Circle-Map Realign`](https://github.com/iprada/Circle-Map)
   1. [`Circle-Map ReadExtractor`](https://github.com/iprada/Circle-Map) -> [`Circle-Map Repeats`](https://github.com/iprada/Circle-Map)
   1. [`CIRCexplorer2`](https://circexplorer2.readthedocs.io/en/latest/)
   1. [`Samblaster`](https://github.com/GregoryFaust/samblaster) -> [`Circle_finder`](https://github.com/pk7zuva/Circle_finder)
   1. Identification of circular amplicons [`AmpliconArchitect`](https://github.com/virajbdeshpande/AmpliconArchitect)
   1. DeNovo Assembly of circular DNAs [`Unicycler`](https://github.com/rrwick/Unicycler) -> [`Minimap2`](https://github.com/lh3/minimap2)
7. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Functionality Overview

A graphical view of the pipeline and its diverse branches can be seen below.

<p align="center">
<img src="docs/images/circdna_pipeline_metromap.png" alt="nf-core/circdna metromap" width="70%">
</p>

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```console
   nextflow run nf-core/circdna -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   ```console
   nextflow run nf-core/circdna --input samplesheet.csv --outdir <OUTDIR> --genome GRCh38 -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
   ```

## Available circular DNA identifiers

Please specify the parameter `circle_identifier` depending on the pipeline branch used for ecDNA identifaction. Please note that some branches/software are only tested with specific NGS data sets.

### Identification of putative ecDNA junctions with ATAC-seq or Circle-seq data

> `circle_finder` uses [Circle_finder](https://github.com/pk7zuva/Circle_finder)

> `circexplorer2` uses [CIRCexplorer2](https://circexplorer2.readthedocs.io/en/latest/)

> `circle_map_realign` uses [Circle-Map Realign](https://github.com/iprada/Circle-Map)

> `circle_map_repeats` uses [Circle-Map Repeats](https://github.com/iprada/Circle-Map) for the identification of repetetive circular DNA

### Identification of circular amplicons with WGS data

> `ampliconarchitect` uses [AmpliconArchitect](https://github.com/virajbdeshpande/AmpliconArchitect)

### De novo assembly of circular DNAs with Circle-Seq data

> `unicycler` uses [Unicycler](https://github.com/rrwick/Unicycler) for de novo assembly of circular DNAs and [Minimap2](https://github.com/lh3/minimap2) for accurate mapping of the identified circular sequences.

### Example Usage

The user can specify either one or multiple `circle_identifier` (see below).

```console
nextflow run nf-core/circdna --input samplesheet.csv --outdir <OUTDIR> --genome GRCh38 -profile docker --circle_identifier circle_map_realign,unicycler
```

## Documentation

The nf-core/circdna pipeline comes with documentation about the pipeline [usage](https://nf-co.re/circdna/usage), [parameters](https://nf-co.re/circdna/parameters) and [output](https://nf-co.re/circdna/output).

## Credits

Main authors:

- [Daniel Schreyer](https://github.com/DSchreyer), University of Glasgow, Institute of Cancer Sciences, Peter Bailey Lab

### Funding

Daniel Schreyer received funding from the European Union’s Horizon 2020 Research and Innovation Program under the Marie Skłodowska-Curie grant agreement No 861196 designated for PRECODE.

nf-core/circdna was originally written by Daniel Schreyer.

<!-- We thank the following people for their extensive assistance in the development of this pipeline: -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#circdna` channel](https://nfcore.slack.com/channels/circdna) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/circdna for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
