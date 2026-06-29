# nf-core/circdna: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2dev - [2026-06-29]

### Enhancements & fixes

- Merged nf-core template updates from v3.0.0 through v4.0.2
  - Updated all CI workflows (linting, nf-test, download pipeline) to nf-core v4 standards
  - Updated GitHub Actions action versions
  - Added ARM64/AMD64 container config files
  - Updated nf-test version to 0.9.4 and Nextflow minimum version requirement
- Migrated from local modules to nf-core modules where available
  - Replaced local BWA index/mem, samtools, trimgalore, fastqc, cat/fastq, unicycler, minimap2 modules with nf-core equivalents
  - Replaced local circexplorer2/parse module with nf-core module
  - Replaced local picard/markduplicates with nf-core `bam_markduplicates_picard` subworkflow
  - Replaced local samtools stats/flagstat/idxstats with nf-core `bam_stats_samtools` subworkflow
- Removed CNVkit local module — CNVkit is already bundled inside AmpliconSuite-Pipeline
- Removed unused local MultiQC module — pipeline uses the nf-core multiqc module
- Removed local CIRCexplorer2 module in favour of nf-core equivalent
- Updated nf-core FastQC module to latest version
- Updated nf-core MultiQC module to latest version (v1.34, using new tuple-based input signature)
- Updated local AmpliconSuite module to use PrepareAA container v1.0.5 (latest)
- Fixed duplicate `$args` bug in the AmpliconSuite-Pipeline script call
- Replaced `Channel` (deprecated) with `channel` factory throughout workflow
- Updated pipeline to use `nf-schema` plugin instead of `nf-validation`
- Updated `utils_nfcore_circdna_pipeline` subworkflow to nf-core v4 structure
  - UTILS_NFSCHEMA_PLUGIN now uses variable-based parameters for help text
- Removed max memory/time/cpu overrides from process configs (use nf-core defaults)
- Removed unused `imNotification` process
- Added nf-test pipeline test (`tests/default.nf.test`) for CI
- Added test dataset (`test_dataset/`) for local development and testing

### Dependencies

| Tool | Previous version | New version |
|------|-----------------|-------------|
| FastQC | 0.11.9 | 0.12.1 |
| MultiQC | 1.18 | 1.34 |
| PrepareAA (AmpliconSuite-Pipeline) | 1.0.3 | 1.0.5 |

## v1.1 - [2024-02-03]

### Credits

Special thanks to the following for their input and contributions to the release:

- [Jens Luebeck](https://github.com/jluebeck)
- [Simon Pearce](https://github.com/SPPearce)
- [Maxime U Garcia](https://github.com/maxulysse)
- [Alex M. Ascensión](https://github.com/alexmascension)

### Enhancements & fixes

- Nf-core template update to 2.11.1
  - update of nf-core modules versions
- Removed AmpliconArchitect and AmpliconClassifier modules with their respective scripts in /bin
  - AmpliconArchitect and AmpliconClassifier is now run inside the AmpliconSuite-Pipeline. Additional scripts are not necessary.
  - Removed respective configs and workflow code
- Added AmpliconSuite-Pipeline
  - A wrapper for calling copy numbers, preparing amplified intervals, running AmpliconArchitect, and calling amplicon classes using AmpliconClassifier
  - Added docker container named [PrepareAA](https://quay.iorepository/nf-core/prepareaa?tab=tags) to run AmpliconSuite-Pipeline with singualarity or docker
  - Added module configs and description
- Changed `assets/multiqc_config.yml`to fit new pipeline version
- Included directory checks for `mosek_license_dir` and `aa_data_repo` .
  - Removed both directory parameters in the test profile as it is only checked when running `ampliconarchitect`
- Updated `nextflow_schema.json` to give better details about how to use `--circle_identifier`
- made `--circle_identifier` an essential parameter
- made `--input_format` an essential parameter and removed the default value to request specification by user
- Updated `--bwa_index` to accept only directory paths to the bwa index files. Makes the user input easier to not need to deal with file endings and patterns. Bug identified by [Alex M. Ascensión](https://github.com/alexmascension) in <https://github.com/nf-core/circdna/issues/68>

## v1.0.4 - [2023-06-26]

### `Added`

### `Fixed`

- Bug that the pipeline only runs with one sample when Picard Markduplicates is used

### `Dependencies`

### `Deprecated`

## v1.0.3 - [2023-05-26]

### `Added`

- Licence, contact, source information for AmpliconArchitect and PrepareAA python scripts
- documentation about absolute path needed of AmpliconArchitect data repository
- ampliconclassifier stub run tests
- new version of circdna metromap with updated colors
- note that ATAC-seq should be used in caution with the pipeline.
- build docker container for prepareaa -> Needs to be built first and will be included in the next release
- nf-core template update 2.8

### `Fixed`

- Circle_finder bug with bash sort command wanting to write into /tmp/ directory and not into work directory
- Usage.md updated to new paths and addition of nf-core modules

### `Dependencies`

### `Deprecated`

- Local python scripts not included in the pipeline
- Local versions of nf-core modules

## v1.0.2 - [2023-03-07]

### `Added`

- ampliconclassifier/makeinput module added -> Generates the input file used for ampliconclassifier functions
- ampliconclassifier/makeresultstable added -> Generates results table from AmpliconArchitect and AmpliconClassifier
- CNN Reference File For AmpliconArchitect
- mm10 option for AmpliconArchitect
- stub runs for AmpliconArchitect processes
- New module versions
- nf-core template 2.7.2

### `Fixed`

- Fixed ZeroDivisionError by Circle-Map
- Fixed keep_duplicates and skip_markduplicates parameter bug

### `Dependencies`

### `Deprecated`

- AmpliconArchitect Summary Process was deprecated

## v1.0.1 - [2022-06-22]

### `Added`

- Documentation Updates

### `Fixed`

- Fixed Bug with pipeline version in nextflow.config
- Fixed Circle-Map Realign bug in which only one sample is processed

### `Dependencies`

### `Deprecated`

- Samtools FAIDX

## v1.0.0 - [2022-06-01]

Initial release of nf-core/circdna, created with the [nf-core](https://nf-co.re/) template.

nf-core/circdna is a bioinformatics analysis pipeline for the identification of circular DNAs in eukaryotic cells. The pipeline is able to process WGS, ATAC-seq data or Circle-Seq data to give insights into the circular DNA landscape in your samples.

In total, the user can choose between 5 different branches inside the pipeline, depending on the biological question and the input data set. In these branches, specific software is used that is built for either the identification of amplified circular DNAs, the detection of putative circular DNA junctions, or the de novo assembly and mapping of circular DNAs.
