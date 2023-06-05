# nf-core/circdna: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [FastQC](#fastqc) - Raw read QC
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution
- [TrimGalore](#trimgalore) - Read Trimming
- [BWA](#bwa) - Read mapping to reference genome
- [Samtools](#samtools) - Sorting, indexing, filtering & stats generation of BAM file
- [Circle-Map Realign](#circle-map-realign) - Identifies putative circular DNA junctions
- [Circle-Map Repeats](#circle-map-repeats) - Identifies putative repetitive circular DNA
- [CIRCexplorer2](#circexplorer2) - Identifies putative circular DNA junctions
- [Circle_finder](#circle_finder) - Identifies putative circular DNA junctions
- [AmpliconArchitect](#ampliconarchitect) - Reconstruct the structure of focally amplified regions
- [Unicycler](#unicycler) - DeNovo Alignment of circular DNAs

## General Tools

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

### TrimGalore

<details markdown="1">
<summary>Output files</summary>

- `trimgalore/`
  - `*_trimming_report.txt`: Trimgalore trimming report.
  - `fastqc/*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.
  - `fastqc/*_fastqc.html`: FastQC report containing quality metrics.

</details>

[TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) combines the trimming tool [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) for the removal of adapter sequences, primers and other unwanted sequences with the quality control tool [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

### BWA

[BWA](https://github.com/lh3/bwa) is a software package for mapping low-divergent sequences against a large reference genome.

Such files are intermediate and not kept in the final files delivered to users.

<details markdown="1">
<summary>Output files</summary>

**Output directory: `results/Reports/[SAMPLE]/SamToolsStats`**

- `[SAMPLE].bam`
  - Alignment file containing information about the read alignment to the reference genome

</details>

### Samtools

#### samtools stats

[samtools stats](https://www.htslib.org/doc/samtools.html) collects statistics from `BAM` files and outputs in a text format.

<details markdown="1">
<summary>Output files</summary>

**Output directory: `results/Reports/[SAMPLE]/SamToolsStats`**

- `[SAMPLE].bam.samtools.stats.out`
  - Raw statistics used by `MultiQC`

Plots will show:

- Alignment metrics.

</details>

For further reading and documentation see the [`samtools` manual](https://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS)

### Mark Duplicates

#### GATK MarkDuplicates

By default, `circdna` will use [GATK MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360042477492-MarkDuplicates-Picard), which locates and tags duplicate reads in a `BAM` or `SAM` file, where duplicate reads are defined as originating from a single fragment of DNA.

<details markdown="1">
<summary>Output files</summary>

**Output directory: `results/markduplicates/bam`**

- `[SAMPLE].md.bam` and `[SAMPLE].md.bai`
  - `BAM` file and index

</details>

For further reading and documentation see the [data pre-processing for variant discovery from the GATK best practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery).

#### Samtools view - Duplicates Filtering

By default, `circdna` removes all duplicates marked by [GATK MarkDuplicates](https://gatk.broadinstitute.org/hc/en-us/articles/360042477492-MarkDuplicates-Picard) using [samtools view](http://www.htslib.org/)

<details markdown="1">
<summary>Output files</summary>

**Output directory: `results/markduplicates/duplicates_removed`**

- `[SAMPLE].md.filtered.sorted.bam` and `[SAMPLE].md.filtered.sorted.bai`
  - `BAM` file and index

</details>

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

## circdna branches

### Branch: `circle_finder`

#### **Circle_finder**

<details markdown="1">
<summary>Output files</summary>

**Output directory: `results/circlefinder/`**

- `[SAMPLE].microDNA-JT.txt`
  - `BED` file containing information about putative circular DNA regions

</details>

[Circle_finder](https://github.com/pk7zuva/Circle_finder) identifies putative circular DNA junctions from paired-end sequencing data. `Circle_finder` uses split and discordant read information to identify junctions that could be generated through the formation of ecDNAs. For more information please see [Circle_finder](https://github.com/pk7zuva/Circle_finder).

### Branch: `circexplorer2`

#### **CIRCexplorer2**

[CIRCexplorer2](https://circexplorer2.readthedocs.io/en/latest/) identifies putative circular DNA junctions from paired-end sequencing data. `CIRCexplorer2` was developed to identify circular RNAs from RNA-seq data. However, it can be also used to call putative circular DNAs from genomic data. For more information see [CIRCexplorer2 docs](https://circexplorer2.readthedocs.io/en/latest/)

<details markdown="1">
<summary>Output files</summary>

**Output directory: `results/circexplorer2/`**

- `[SAMPLE].circexplorer_circdna.bed`
  - `BED` file containing information about putative circular DNA regions
- `[SAMPLE].CIRCexplorer2_parse.log`
  - `log` file

</details>

### Branch: `circle_map_realign`

`circle_map_realign` uses the functionality of `Circle-Map` to call putative circular DNAs from mappable regions. To identify circular DNAs it uses information about split and discordant reads and uses realignment steps to identify the exact breakpoint of the circular DNA. For more information, please see the [original paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3160-3) or the [GitHub Page](https://github.com/iprada/Circle-Map)

#### **Circle-Map Readextractor**

[Circle-Map Readextractor](https://github.com/iprada/Circle-Map) extracts read candidates for circular DNA identification.

<details markdown="1">
<summary>Output files</summary>

**Output directory: `results/circlemap/readextractor`**

- `[SAMPLE].qname.sorted.circular_read_candidates.bam`
  - `BAM` file containing candidate reads

</details>

#### **Circle-Map Realign**

[Circle-Map Realign](https://github.com/iprada/Circle-Map) detects putative circular DNA junctions from read candidates extracted by `Circle-Map Readextractor`

<details markdown="1">
<summary>Output files</summary>

**Output directory: `results/circlemap/realign`**

- `[SAMPLE]_circularDNA_coordinates.bed`
  - `BED` file containing information about putative circular DNA regions

</details>

### Branch: `circle_map_repeats`

#### **Circle-Map Readextractor**

[Circle-Map Readextractor](https://github.com/iprada/Circle-Map) extracts read candidates for circular DNA identification.

<details markdown="1">
<summary>Output files</summary>

**Output directory: `results/circlemap/readextractor`**

- `[SAMPLE].qname.sorted.circular_read_candidates.bam`
  - `BAM` file containing candidate reads

</details>

#### **Circle-Map Repeats**

[Circle-Map Repeats](https://github.com/iprada/Circle-Map) identifies chromosomal coordinates from repetetive circular DNAs.

<details markdown="1">
<summary>Output files</summary>

**Output directory: `results/circlemap/repeats`**

- `[SAMPLE]_circularDNA_repeats_coordinates.bed`
  - `BED` file containing information about repetetive circular DNAs

</details>

### Branch: `unicycler`

This Branch utilises the ability of [Unicycler](https://github.com/rrwick/Unicycler) to denovo assemble circular DNAs in combination with the long read mapping capabilities of [Minimap2](https://github.com/lh3/minimap2), to identify the origin of the circular DNAs.

#### **Unicycler**

[Unicycler](https://github.com/rrwick/Unicycler) was originally built as an assembly pipeline for bacterial genomes. In `nf-core/circdna` it is used to denovo assemble circular DNAs.

<details markdown="1">
<summary>Output files</summary>

**Output directory: `results/unicycler/`**

- `[SAMPLE].assembly.gfa.gz`
  - `gfa` file containing sequence of denovo assembled sequences
- `[SAMPLE].assembly.scaffolds.fa.gz`
  - `fasta` file containing sequences of denovo assembled sequences in fasta format with information if denovo assembled seoriginated from a circular DNA.quence forms a circular contig.

</details>

#### **Minimap2**

[Minimap2](https://github.com/lh3/minimap2) uses circular DNA sequences identified by Unicycler and maps it to the given reference genome.

<details markdown="1">
<summary>Output files</summary>

**Output directory: `results/unicycler/minimap2`**

- `[SAMPLE].paf`
  - `paf` file containing mapping information of circular DNA sequences

</details>

### Branch: `ampliconarchitect`

This pipeline branch `ampliconarchitect` is only usable with WGS data. This branch uses the utility of [PrepareAA](https://github.com/jluebeck/Prepare) to collect amplified seeds from copy number calls, which will be then fed to [AmpliconArchitect](https://github.com/jluebeck/AmpliconArchitect) to characterise amplicons in each given sample.

#### **CNVkit**

[CNVkit](https://cnvkit.readthedocs.io/en/stable/) uses alignment information to make copy number calls. These copy number calls will be used by AmpliconArchitect to identify circular and other types of amplicons. The Copy Number calls are then connected to seeds and filtered based on the copy number threshold using scripts provided by [PrepareAA](https://github.com/jluebeck/AmpliconSuite-pipeline)

<details markdown="1">
<summary>Output files</summary>

**Output directory: `results/cnvkit`**

- `[SAMPLE]_CNV_GAIN.bed`
  - `bed` file containing filtered Copy Number calls
- `[SAMPLE]_AA_CNV_SEEDS.bed`
  - `bed` file containing filtered and connected amplified regions (seeds). This is used as input for [AmpliconArchitect](https://github.com/jluebeck/AmpliconArchitect)
- `[SAMPLE].cnvkit.segment.cns`
  - `cns` file containing copy number calls of CNVkit segment.
- `[SAMPLE].cnvkit.segment.cnr`
  - `cns` file containing copy number calls of CNVkit in cnr format.

</details>

#### **AmpliconArchitect**

[AmpliconArchitect](https://github.com/jluebeck/AmpliconArchitect) uses amplicon seeds provided by `CNVkit`and `PrepareAA`to identify different types of amplicons in each sample.

<details markdown="1">
<summary>Output files</summary>

**Output directory: `results/ampliconarchitect/ampliconarchitect`**

- `amplicons/[SAMPLE]_[AMPLICONID]_cycles.txt`
  - `txt`file describing the amplicon segments
- `amplicons/[SAMPLE]_[AMPLICONID]_graph.txt`
  - `txt` file describing the amplicon graph
- `cnseg/[SAMPLE]_[SEGMENT]_graph.txt`
  - `txt` file describing the copy number segmentation file
- `summary/[SAMPLE]_summary.txt`
  - `txt` file describing each amplicon with regards to breakpoints, composition, oncogene content, copy number
- `sv_view/[SAMPLE]_[AMPLICONID].{png,pdf}`
  - `png` or `pdf` file displaying the amplicon rearrangement signature

</details>

#### **AmpliconClassifier**

[AmpliconClassifier](https://github.com/jluebeck/AmpliconClassifier) classifies each amplicon by using the `cycles` and the `graph` files generated by `AmpliconArchitect`.

<details markdown="1">
<summary>Output files</summary>

**Output directory: `results/ampliconclassifier`**

- `makeinput/ampliconclassifier.input`
  - `txt` file containing the input used for `AmpliconClassifier` and `AmpliconSimilarity`.
- `ampliconclassifier/ampliconclassifier_amplicon_classification_profiles.tsv`
  - `tsv` file describing the amplicon class of each amplicon for each sample.
- `ecDNA_counts/ampliconclassifier_ecDNA_counts.tsv`
  - `tsv` file describing if an amplicon is circular [1 = circular, 0 = non-circular].
- `gene_list/ampliconclassifier_gene_list.tsv`
  - `tsv` file detailing the genes on each amplicon.
- `log/ampliconclassifier_stdout.log`
  - `log` file
- `ampliconsimilarity/ampliconclassifier_similarity_scores.tsv`
  - `tsv` file containing amplicon similarity scores calculated by `AmpliconSimilarity`.
- `bed/[SAMPLE]_amplicon[AMPLICONID]_[CLASSIFICATION]_[ID]_intervals.bed`
  - `bed` files containing information about the intervals on each amplicon. `unknown` intervals were not identified to be located on the respective amplicon.
- `resultstable/ampliconclassifier_result_table.[tsv,json]`
  - `tsv` or `json` file of the results table tenerated by `AmpliconClassifier` which combines the output of `AmpliconArchitect` and `AmpliconClassifier`.

</details>

### Reports

#### samtools stats

[samtools stats](https://www.htslib.org/doc/samtools.html) collects statistics from CRAM files and outputs in a text format.
For further reading and documentation see the [`samtools` manual](https://www.htslib.org/doc/samtools.html#COMMANDS_AND_OPTIONS).

The plots will show:

- Alignment metrics.

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/reports/samtools_stats/<sample>`**

- `[SAMPLE].flagstat`
  - Samtools flagstat outpu
- `[SAMPLE].idxstats`
  - Samtools idxstats output
- `[SAMPLE].stats`
  - Samtools stats output

</details>
### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
