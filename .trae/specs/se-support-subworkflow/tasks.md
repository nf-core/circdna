# SE Support & Subworkflow Splitting - Implementation Plan

## [x] Task 1: Update samplesheet parsing for SE support

- **Priority**: high
- **Depends On**: None
- **Description**:
  - Modify `subworkflows/local/input_check` to handle samplesheets without `fastq_2` column
  - Add logic to detect SE vs PE data
  - Update `nextflow_schema.json` to make `fastq_2` optional
- **Acceptance Criteria Addressed**: AC-1, AC-2
- **Test Requirements**:
  - `programmatic` TR-1.1: Parse SE samplesheet and create channel with single fastq file
  - `programmatic` TR-1.2: Parse PE samplesheet and create channel with two fastq files
- **Notes**: Need to check how nf-validation plugin handles optional columns

## [x] Task 2: Update CAT_FASTQ for SE support

- **Priority**: high
- **Depends On**: Task 1
- **Description**:
  - Modify CAT_FASTQ module usage in main workflow to handle SE data
  - Ensure proper channel handling for single-file inputs
- **Acceptance Criteria Addressed**: AC-3
- **Test Requirements**:
  - `programmatic` TR-2.1: CAT_FASTQ processes SE data without errors
  - `programmatic` TR-2.2: CAT_FASTQ continues to work with PE data
- **Notes**: CAT_FASTQ module should already support single-file inputs

## [x] Task 3: Update TRIMGALORE for SE support

- **Priority**: high
- **Depends On**: Task 1
- **Description**:
  - Modify TRIMGALORE module call to add `--single` flag for SE data
  - Handle single fastq file input properly
- **Acceptance Criteria Addressed**: AC-3
- **Test Requirements**:
  - `programmatic` TR-3.1: TRIMGALORE processes SE data with --single flag
  - `programmatic` TR-3.2: TRIMGALORE continues to work with PE data
- **Notes**: Need to check nf-core TRIMGALORE module interface

## [x] Task 4: Update BWA_MEM for SE support

- **Priority**: high
- **Depends On**: Task 1
- **Description**:
  - Modify BWA_MEM module call to handle SE data
  - BWA_MEM needs to know whether data is SE or PE
- **Acceptance Criteria Addressed**: AC-3
- **Test Requirements**:
  - `programmatic` TR-4.1: BWA_MEM aligns SE data correctly
  - `programmatic` TR-4.2: BWA_MEM continues to work with PE data
- **Notes**: BWA_MEM module has a `paired` parameter that needs to be set

## [x] Task 5: Update CIRCLEMAP_READEXTRACTOR for SE support

- **Priority**: high
- **Depends On**: Task 1
- **Description**:
  - Check and update local CIRCLEMAP_READEXTRACTOR module for SE support
  - Ensure it can handle SE BAM files
- **Acceptance Criteria Addressed**: AC-3
- **Test Requirements**:
  - `programmatic` TR-5.1: CIRCLEMAP_READEXTRACTOR processes SE data without errors
- **Notes**: This is a local module, may need code changes

## [x] Task 6: Create bam_preprocessing subworkflow

- **Priority**: medium
- **Depends On**: None
- **Description**:
  - Create `subworkflows/local/bam_preprocessing/main.nf`
  - Include: BWA_MEM, SAMTOOLS_INDEX, BAM_STATS_SAMTOOLS, BAM_MARKDUPLICATES_PICARD
  - Define clear input/output channels
- **Acceptance Criteria Addressed**: AC-4
- **Test Requirements**:
  - `programmatic` TR-6.1: bam_preprocessing subworkflow runs without errors
  - `human-judgment` TR-6.2: Code follows nf-core style guidelines
- **Notes**: Need to handle markduplicates filtering options

## [x] Task 7: Create circle_map_pipeline subworkflow

- **Priority**: medium
- **Depends On**: Task 6
- **Description**:
  - Create `subworkflows/local/circle_map_pipeline/main.nf`
  - Include: SAMTOOLS_SORT_QNAME, CIRCLEMAP_READEXTRACTOR, SAMTOOLS_SORT_RE, CIRCLEMAP_REALIGN, CIRCLEMAP_REPEATS
- **Acceptance Criteria Addressed**: AC-4
- **Test Requirements**:
  - `programmatic` TR-7.1: circle_map_pipeline subworkflow runs without errors
  - `human-judgment` TR-7.2: Code follows nf-core style guidelines
- **Notes**: Need to handle both realign and repeats branches

## [x] Task 8: Create circle_finder_pipeline subworkflow

- **Priority**: medium
- **Depends On**: Task 6
- **Description**:
  - Create `subworkflows/local/circle_finder_pipeline/main.nf`
  - Include: SAMTOOLS_SORT_QNAME, SAMBLASTER, BEDTOOLS_SPLITBAM2BED, BEDTOOLS_SORTEDBAM2BED, CIRCLEFINDER
- **Acceptance Criteria Addressed**: AC-4
- **Test Requirements**:
  - `programmatic` TR-8.1: circle_finder_pipeline subworkflow runs without errors
  - `human-judgment` TR-8.2: Code follows nf-core style guidelines
- **Notes**: Need to handle multiple BEDTOOLS outputs

## [x] Task 9: Create ampliconarchitect_pipeline subworkflow

- **Priority**: medium
- **Depends On**: Task 6
- **Description**:
  - Create `subworkflows/local/ampliconarchitect_pipeline/main.nf`
  - Include: CNVKIT_BATCH, CNVKIT_SEGMENT, AMPLICONSUITE
- **Acceptance Criteria Addressed**: AC-4
- **Test Requirements**:
  - `programmatic` TR-9.1: ampliconarchitect_pipeline subworkflow runs without errors
  - `human-judgment` TR-9.2: Code follows nf-core style guidelines
- **Notes**: Need to handle CNVkit reference input

## [x] Task 10: Create unicycler_pipeline subworkflow

- **Priority**: low
- **Depends On**: None
- **Description**:
  - Create `subworkflows/local/unicycler_pipeline/main.nf`
  - Include: UNICYCLER, SEQTK_SEQ, GETCIRCULARREADS, MINIMAP2_ALIGN
- **Acceptance Criteria Addressed**: AC-4
- **Test Requirements**:
  - `programmatic` TR-10.1: unicycler_pipeline subworkflow runs without errors
  - `human-judgment` TR-10.2: Code follows nf-core style guidelines
- **Notes**: Unicycler already handles SE/PE automatically

## [x] Task 11: Update main workflow to use subworkflows

- **Priority**: high
- **Depends On**: Tasks 6-10
- **Description**:
  - Modify `workflows/circdna.nf` to call subworkflows instead of individual modules
  - Update channel connections between subworkflows
- **Acceptance Criteria Addressed**: AC-4, AC-5
- **Test Requirements**:
  - `programmatic` TR-11.1: Pipeline runs successfully with subworkflows
  - `human-judgment` TR-11.2: Main workflow code is cleaner and more readable
- **Notes**: This is the most critical task - ensure all channels are properly connected

## [x] Task 12: Update CHANGELOG and CHANGES&FIX documentation

- **Priority**: medium
- **Depends On**: All tasks
- **Description**:
  - Update `CHANGELOG.md` with SE support and subworkflow splitting changes
  - Create new CHANGES&FIX document for this release
- **Acceptance Criteria Addressed**: Documentation requirement
- **Test Requirements**:
  - `human-judgment` TR-12.1: CHANGELOG is updated with all changes
  - `human-judgment` TR-12.2: CHANGES&FIX document is complete
