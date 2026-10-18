# SE Support & Subworkflow Splitting - Product Requirement Document

## Overview

- **Summary**: Add Single-End (SE) sequencing data support to the nf-core/circdna pipeline and split the monolithic main workflow into modular subworkflows for better maintainability and reusability.
- **Purpose**: Enable processing of single-end sequencing data (currently only paired-end is supported) and improve code organization through subworkflow decomposition.
- **Target Users**: Bioinformaticians using single-end sequencing data for circular DNA analysis and pipeline developers maintaining the nf-core/circdna pipeline.

## Goals

- Enable SE data processing throughout the pipeline
- Split monolithic workflow into modular subworkflows
- Maintain backward compatibility with existing PE data processing
- Improve code maintainability and readability

## Non-Goals (Out of Scope)

- Add new circular DNA detection tools
- Modify existing module logic (only adapt for SE support)
- Add new pipeline branches
- Update documentation beyond what's necessary for implementation

## Background & Context

- Current pipeline only supports paired-end (PE) sequencing data
- Main workflow file (`workflows/circdna.nf`) contains ~500 lines of code with all module calls in one place
- SE support is a common requirement for many sequencing experiments
- Subworkflow splitting follows nf-core best practices for code organization

## Functional Requirements

- **FR-1**: Support SE data in samplesheet parsing
- **FR-2**: Support SE data in CAT_FASTQ module
- **FR-3**: Support SE data in TRIMGALORE module
- **FR-4**: Support SE data in BWA_MEM alignment
- **FR-5**: Support SE data in CIRCLEMAP_READEXTRACTOR
- **FR-6**: Split BAM preprocessing into `bam_preprocessing` subworkflow
- **FR-7**: Split Circle-Map pipeline into `circle_map_pipeline` subworkflow
- **FR-8**: Split Circle-Finder pipeline into `circle_finder_pipeline` subworkflow
- **FR-9**: Split AmpliconArchitect pipeline into `ampliconarchitect_pipeline` subworkflow
- **FR-10**: Split Unicycler pipeline into `unicycler_pipeline` subworkflow

## Non-Functional Requirements

- **NFR-1**: Maintain backward compatibility - PE data processing should continue to work
- **NFR-2**: Subworkflows should have clear input/output contracts
- **NFR-3**: Code should follow nf-core style guidelines

## Constraints

- **Technical**: Must use Nextflow DSL2 syntax
- **Dependencies**: Must maintain compatibility with existing nf-core modules
- **Testing**: Must pass existing tests

## Assumptions

- Users will provide properly formatted samplesheets with or without `fastq_2` column
- SE data follows standard FASTQ format
- Existing modules support SE mode with minor configuration changes

## Acceptance Criteria

### AC-1: SE Samplesheet Parsing

- **Given**: A samplesheet with only `sample` and `fastq_1` columns
- **When**: The pipeline is run with this samplesheet
- **Then**: The pipeline should parse the samplesheet correctly and identify it as SE data
- **Verification**: `programmatic`

### AC-2: PE Samplesheet Parsing (Backward Compatibility)

- **Given**: A samplesheet with `sample`, `fastq_1`, and `fastq_2` columns
- **When**: The pipeline is run with this samplesheet
- **Then**: The pipeline should parse the samplesheet correctly and identify it as PE data
- **Verification**: `programmatic`

### AC-3: SE Data Processing

- **Given**: SE sequencing data
- **When**: The pipeline processes the data through all modules
- **Then**: All modules should complete successfully without errors
- **Verification**: `programmatic`

### AC-4: Subworkflow Structure

- **Given**: The main workflow file
- **When**: Examining the code structure
- **Then**: The main workflow should call subworkflows instead of individual modules
- **Verification**: `human-judgment`

### AC-5: Output Compatibility

- **Given**: Both SE and PE data
- **When**: Running the pipeline with both data types
- **Then**: The output files should have the same structure and naming conventions
- **Verification**: `human-judgment`

## Open Questions

- [ ] How should the samplesheet schema be updated to support optional `fastq_2`?
- [ ] Are there any modules that fundamentally don't support SE data?
