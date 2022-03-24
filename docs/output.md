# nf-core/circdna: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

<<<<<<< HEAD
=======
<<<<<<< HEAD
<<<<<<< HEAD
* [FastQC](#fastqc) - Raw read QC
* [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
* [TrimGalore](#trimgalore) - Read Trimming
* [BWA](#bwa) - Read mapping to reference genome
* [Samtools](#samtools) - Sorting, indexing, filtering & stats generation of BAM file
* [Bedtools](#bedtools) - Converts bam to bed file
* [Circle-Map](#circle-map) - Identifies putative circular DNA junctions
* [CIRCexplorer2](#circexplorer2) - Identifies putative circular DNA junctions
* [Circle_finder](#circle_finder) - Identifies putative circular DNA junctions
* [AmpliconArchitect](#ampliconarchitect) - Reconstruct the structure of focally amplified regions
* [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution
=======
<<<<<<< HEAD
-   [FastQC](#fastqc) - Raw read QC
-   [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
-   [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution
=======
>>>>>>> 6a32be9b78b05990be38e99304c76909cdbea507
- [FastQC](#fastqc) - Raw read QC
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [TrimGalore](#trimgalore) - Read Trimming
- [BWA](#bwa) - Read mapping to reference genome
- [Samtools](#samtools) - Sorting, indexing, filtering & stats generation of BAM file
- [Bedtools](#bedtools) - Converts bam to bed file
- [Circle-Map](#circle-map) - Identifies putative circular DNA junctions
- [CIRCexplorer2](#circexplorer2) - Identifies putative circular DNA junctions
- [Circle_finder](#circle_finder) - Identifies putative circular DNA junctions
- [AmpliconArchitect](#ampliconarchitect) - Reconstruct the structure of focally amplified regions
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution
<<<<<<< HEAD
=======
>>>>>>> d13c279908de1b8cc2914a29996b39dc584e9e3f
>>>>>>> nf-core-TEMPLATE
=======
-   [FastQC](#fastqc) - Raw read QC
-   [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
-   [TrimGalore](#trimgalore) - Read Trimming
-   [BWA](#bwa) - Read mapping to reference genome
-   [Samtools](#samtools) - Sorting, indexing, filtering & stats generation of BAM file
-   [Bedtools](#bedtools) - Converts bam to bed file
-   [Circle-Map](#circle-map) - Identifies putative circular DNA junctions
-   [CIRCexplorer2](#circexplorer2) - Identifies putative circular DNA junctions
-   [Circle_finder](#circle_finder) - Identifies putative circular DNA junctions
-   [AmpliconArchitect](#ampliconarchitect) - Reconstruct the structure of focally amplified regions
-   [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution
>>>>>>> 417c56964a4eef354736058b448d677d04172201
>>>>>>> 6a32be9b78b05990be38e99304c76909cdbea507

### FastQC

<details markdown="1">
<summary>Output files</summary>

<<<<<<< HEAD
- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.
=======
<<<<<<< HEAD
<<<<<<< HEAD
-   `fastqc/`
    -   `*_fastqc.html`: FastQC report containing quality metrics.
    -   `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.
=======
- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.
>>>>>>> d13c279908de1b8cc2914a29996b39dc584e9e3f
=======
-   `fastqc/`
    -   `*_fastqc.html`: FastQC report containing quality metrics.
    -   `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.
>>>>>>> 417c56964a4eef354736058b448d677d04172201
>>>>>>> 6a32be9b78b05990be38e99304c76909cdbea507

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

<<<<<<< HEAD
=======
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 417c56964a4eef354736058b448d677d04172201
-   `multiqc/`
    -   `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
    -   `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
    -   `multiqc_plots/`: directory containing static images from the report in various formats.
<<<<<<< HEAD
=======
>>>>>>> 6a32be9b78b05990be38e99304c76909cdbea507
- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.
<<<<<<< HEAD
=======
>>>>>>> d13c279908de1b8cc2914a29996b39dc584e9e3f
=======
>>>>>>> 417c56964a4eef354736058b448d677d04172201
>>>>>>> 6a32be9b78b05990be38e99304c76909cdbea507

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

<<<<<<< HEAD
=======
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 417c56964a4eef354736058b448d677d04172201
-   `pipeline_info/`
    -   Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
    -   Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
    -   Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
<<<<<<< HEAD
=======
>>>>>>> 6a32be9b78b05990be38e99304c76909cdbea507
- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
<<<<<<< HEAD
=======
>>>>>>> d13c279908de1b8cc2914a29996b39dc584e9e3f
=======
>>>>>>> 417c56964a4eef354736058b448d677d04172201
>>>>>>> 6a32be9b78b05990be38e99304c76909cdbea507

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
